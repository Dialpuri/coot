#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"
#include "api/coot-molecule.hh"
#include "lidia-core/rdkit-interface.hh"
#include "geometry/protein-geometry.hh"
#include <rdkit/GraphMol/RDKitBase.h>

static std::string chiral_type_str(RDKit::Atom::ChiralType ct) {
    if (ct == RDKit::Atom::CHI_UNSPECIFIED) return "CHI_UNSPECIFIED";
    if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW) return "CHI_TETRAHEDRAL_CW";
    if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW) return "CHI_TETRAHEDRAL_CCW";
    if (ct == RDKit::Atom::CHI_OTHER) return "CHI_OTHER";
    return "UNKNOWN";
}

// Pad atom name to 4 chars (MMDB/PDB-style)
static std::string pad_name_4c(const std::string& name) {
    std::string s = name;
    if (s.size() >= 4) return s.substr(0, 4);
    // PDB convention: if first char is alpha and length is 1 or 2, left-pad by 1
    if (!s.empty() && std::isalpha(static_cast<unsigned char>(s[0])) && s.size() <= 2) {
        s = " " + s;
    }
    // Right-pad to 4
    while (s.size() < 4) s += ' ';
    return s;
}

TEST(OracleTest, set_atom_chirality) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_GE(imol, 0);

    coot::protein_geometry &geom = mc.get_geometry();
    int imol_enc = mc.get_imol_enc_any();

    // Load the same structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Get restraints for THR (chiral amino acid with known chiral center at C beta)
    {
        std::pair<bool, coot::dictionary_residue_restraints_t> restraints_pair =
            geom.get_monomer_restraints("THR", imol_enc);

        EXPECT_TRUE(restraints_pair.first);
        if (restraints_pair.first) {
            const coot::dictionary_residue_restraints_t &restraints = restraints_pair.second;
            EXPECT_EQ(restraints.chiral_restraint.size(), 2);
            EXPECT_EQ(restraints.atom_info.size(), 17);

            // Find gemmi::Residue and gemmi::Atom for chain A, residue 21
            const gemmi::Residue* res_ptr = nullptr;
            const gemmi::Atom* at_ptr = nullptr;

            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        if (res.seqid.num.value == 21) {
                            res_ptr = &res;
                            // Pick first non-H atom
                            for (gemmi::Atom& atom : res.atoms) {
                                if (!atom.element.is_hydrogen()) {
                                    at_ptr = &atom;
                                    break;
                                }
                            }
                            if (!at_ptr && !res.atoms.empty()) {
                                at_ptr = &res.atoms[0];
                            }
                            break;
                        }
                    }
                    break;
                }
            }
            ASSERT_NE(res_ptr, nullptr);
            ASSERT_NE(at_ptr, nullptr);

            // Case 1a: atom whose name matches a chiral restraint center
            if (restraints.chiral_restraint.size() > 0) {
                std::string chiral_atom_name = restraints.chiral_restraint[0].atom_id_c_4c();
                EXPECT_EQ(chiral_atom_name, " CA ");

                // Find matching gemmi::Atom with this name (strip whitespace for comparison)
                std::string chiral_atom_name_stripped = coot::util::remove_whitespace(chiral_atom_name);
                for (const gemmi::Atom& atom : res_ptr->atoms) {
                    std::string atom_name_stripped = coot::util::remove_whitespace(atom.name);
                    if (atom_name_stripped == chiral_atom_name_stripped) {
                        at_ptr = &atom;
                        break;
                    }
                }

                RDKit::Atom *rdkit_at = new RDKit::Atom;
                const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
                // Use gemmi element name
                std::string ele = at_ptr->element.name();
                int z = tbl->getAtomicNumber(ele);
                rdkit_at->setAtomicNum(z);

                RDKit::Atom::ChiralType before = rdkit_at->getChiralTag();
                EXPECT_EQ(chiral_type_str(before), "CHI_UNSPECIFIED");

                coot::set_atom_chirality_gemmi(rdkit_at, at_ptr, res_ptr, restraints);

                RDKit::Atom::ChiralType after = rdkit_at->getChiralTag();
                EXPECT_EQ(chiral_type_str(after), "CHI_TETRAHEDRAL_CCW");

                // Check properties
                std::string prop_val;
                try {
                    rdkit_at->getProp("mmcif_chiral_N1", prop_val);
                    EXPECT_EQ(prop_val, "N");
                } catch (...) {
                    FAIL() << "mmcif_chiral_N1 not set";
                }
                try {
                    rdkit_at->getProp("mmcif_chiral_volume_sign", prop_val);
                    EXPECT_EQ(prop_val, "positive");
                } catch (...) {
                    FAIL() << "mmcif_chiral_volume_sign not set";
                }

                delete rdkit_at;
            }

            // Case 1b: atom NOT a chiral center (but still CA in this case)
            {
                RDKit::Atom *rdkit_at = new RDKit::Atom;
                const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
                std::string ele = at_ptr->element.name();
                int z = tbl->getAtomicNumber(ele);
                rdkit_at->setAtomicNum(z);

                EXPECT_EQ(pad_name_4c(at_ptr->name), " CA ");

                RDKit::Atom::ChiralType before = rdkit_at->getChiralTag();
                EXPECT_EQ(chiral_type_str(before), "CHI_UNSPECIFIED");

                coot::set_atom_chirality_gemmi(rdkit_at, at_ptr, res_ptr, restraints);

                RDKit::Atom::ChiralType after = rdkit_at->getChiralTag();
                EXPECT_EQ(chiral_type_str(after), "CHI_TETRAHEDRAL_CCW");

                delete rdkit_at;
            }
        }
    }

    // Case 2: empty restraints (no chiral info at all)
    {
        coot::dictionary_residue_restraints_t empty_restraints;

        // Find chain A, residue 21, first atom from gemmi
        const gemmi::Atom* at_ptr = nullptr;
        const gemmi::Residue* res_ptr = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 21) {
                        res_ptr = &res;
                        if (!res.atoms.empty()) {
                            at_ptr = &res.atoms[0];
                        }
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(at_ptr, nullptr);

        RDKit::Atom *rdkit_at = new RDKit::Atom;
        const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
        rdkit_at->setAtomicNum(6); // carbon

        EXPECT_EQ(chiral_type_str(rdkit_at->getChiralTag()), "CHI_UNSPECIFIED");

        coot::set_atom_chirality_gemmi(rdkit_at, at_ptr, res_ptr, empty_restraints);

        EXPECT_EQ(chiral_type_str(rdkit_at->getChiralTag()), "CHI_UNSPECIFIED");

        delete rdkit_at;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
