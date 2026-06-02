#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, water_coordination_t_add_contact_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u) << "Expected at least one model";

    gemmi::Model& model = st.models[0];

    // Get atoms: water O from chain B, and CA from chain A
    gemmi::Atom const* water_O = nullptr;
    gemmi::Atom const* protein_CA = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                if (chain.name == "B" && res.seqid.num.value == 1 && atom.name == "O") {
                    water_O = &atom;
                }
                if (chain.name == "A" && res.seqid.num.value == 1 && atom.name == "CA") {
                    protein_CA = &atom;
                }
            }
        }
    }
    ASSERT_NE(water_O, nullptr) << "water_O not found";
    ASSERT_NE(protein_CA, nullptr) << "protein_CA not found";

    // Identity matrix
    double mat[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            mat[i][j] = (i == j) ? 1.0 : 0.0;

    // Case 1: valid contact (both have empty altLoc -> ' ', protein element != "H")
    {
        coot::util::water_coordination_t wc;
        size_t before = wc.atom_contacts.size();
        wc.add_contact_gemmi(water_O, protein_CA, mat);
        size_t after = wc.atom_contacts.size();

        EXPECT_EQ(before, 0u);
        EXPECT_EQ(after, 1u);
    }

    // Case 2: hydrogen contactor should be filtered out
    // Find a hydrogen atom in the structure to use as contactor
    {
        gemmi::Atom const* hydrogen_atom = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (atom.element.is_hydrogen()) {
                        hydrogen_atom = &atom;
                        break;
                    }
                }
                if (hydrogen_atom) break;
            }
            if (hydrogen_atom) break;
        }

        if (hydrogen_atom) {
            coot::util::water_coordination_t wc2;
            size_t before2 = wc2.atom_contacts.size();
            wc2.add_contact_gemmi(water_O, hydrogen_atom, mat);
            size_t after2 = wc2.atom_contacts.size();

            EXPECT_EQ(before2, 0u);
            EXPECT_EQ(after2, 0u);
        }
    }

    // Case 3: non-matching altLocs should be rejected
    // Need atoms with non-matching non-space altlocs.
    // We can't modify gemmi atoms, so look for existing ones with altlocs.
    {
        std::vector<gemmi::Atom const*> altA_atoms;
        std::vector<gemmi::Atom const*> altB_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (atom.altloc == 'A') altA_atoms.push_back(&atom);
                    if (atom.altloc == 'B') altB_atoms.push_back(&atom);
                }
            }
        }

        if (!altA_atoms.empty() && !altB_atoms.empty()) {
            coot::util::water_coordination_t wc3;
            size_t before3 = wc3.atom_contacts.size();
            wc3.add_contact_gemmi(altA_atoms[0], altB_atoms[0], mat);
            size_t after3 = wc3.atom_contacts.size();

            EXPECT_EQ(before3, 0u);
            EXPECT_EQ(after3, 0u);
        }
    }

    // Case 4: matching non-empty altLocs should produce a contact
    {
        std::vector<gemmi::Atom const*> altA_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (atom.altloc == 'A') altA_atoms.push_back(&atom);
                }
            }
        }

        if (altA_atoms.size() >= 2) {
            coot::util::water_coordination_t wc4;
            size_t before4 = wc4.atom_contacts.size();
            wc4.add_contact_gemmi(altA_atoms[0], altA_atoms[1], mat);
            size_t after4 = wc4.atom_contacts.size();

            EXPECT_EQ(before4, 0u);
            EXPECT_EQ(after4, 1u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}