#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_chiral_tag_v2) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];

    // Get residue VAL A 7 (has chiral center at CA)
    gemmi::Residue *res = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue &residue : chain.residues) {
                if (residue.name == "VAL" && residue.seqid.num.value == 7) {
                    res = &residue;
                    break;
                }
            }
        }
    }
    EXPECT_NE(res, nullptr);

    // Get atom CA from VAL
    const gemmi::Atom *atom_CA = nullptr;
    for (const gemmi::Atom &atom : res->atoms) {
        std::string aname = atom.name;
        std::cout << "Atom name: '" << aname << "' len=" << aname.size() << std::endl;
        if (aname == "CA") {
            atom_CA = &atom;
            break;
        }
    }
    EXPECT_NE(atom_CA, nullptr);

    // Case 1: Empty restraints — should return CHI_UNSPECIFIED (0)
    {
        coot::dictionary_residue_restraints_t empty_restraints;
        RDKit::Atom::ChiralType tag = coot::get_chiral_tag_v2_gemmi(res, empty_restraints, atom_CA);
        EXPECT_EQ(static_cast<int>(tag), 0);
    }

    // Case 2: Full chiral + bond restraints for VAL CA
    // CA is chiral centre bonded to N, CB, C (from chiral restraint) + HA (4th via bond)
    {
        coot::dictionary_residue_restraints_t restraints;

        // Chiral restraint: CA=centre, N/CB/C=3 chiral neighbours
        coot::dict_chiral_restraint_t cr("VAL_chiral", "CA", "N", "CB", "C", 1);
        restraints.chiral_restraint.push_back(cr);

        // Bond restraints: CA bonded to N, C, CB, HA
        coot::dict_bond_restraint_t b1, b2, b3, b4;
        b1.set_atom_id_1("CA"); b1.set_atom_id_2("N");
        b2.set_atom_id_1("CA"); b2.set_atom_id_2("C");
        b3.set_atom_id_1("CA"); b3.set_atom_id_2("CB");
        b4.set_atom_id_1("CA"); b4.set_atom_id_2("HA");
        restraints.bond_restraint.push_back(b1);
        restraints.bond_restraint.push_back(b2);
        restraints.bond_restraint.push_back(b3);
        restraints.bond_restraint.push_back(b4);

        RDKit::Atom::ChiralType tag = coot::get_chiral_tag_v2_gemmi(res, restraints, atom_CA);
        EXPECT_EQ(static_cast<int>(tag), 2);
    }

    // Case 3: Non-existent chiral centre — should return CHI_UNSPECIFIED (0)
    {
        coot::dictionary_residue_restraints_t restraints;
        coot::dict_chiral_restraint_t cr("nonexist", "NONEXISTENT", "N", "CB", "C", 1);
        restraints.chiral_restraint.push_back(cr);

        RDKit::Atom::ChiralType tag = coot::get_chiral_tag_v2_gemmi(res, restraints, atom_CA);
        EXPECT_EQ(static_cast<int>(tag), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}