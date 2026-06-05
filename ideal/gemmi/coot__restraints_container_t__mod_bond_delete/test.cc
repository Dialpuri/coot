#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mod_bond_delete) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find residue A/21 (ALA)
    gemmi::Model &model = st.models[0];
    const gemmi::Residue *residue_p = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue &res : chain.residues) {
                if (res.seqid.num.value == 21 && res.name == "ALA") {
                    residue_p = &res;
                    break;
                }
            }
        }
    }
    ASSERT_TRUE(residue_p != nullptr);

    // Build a simplified restraints_container_t with 6 BOND restraints for ALA/21
    coot::restraints_container_t restraints;

    // Populate atoms from the residue using padded_name() for 4-char PDB format
    for (const gemmi::Atom &atom : residue_p->atoms) {
        restraints.atoms.push_back({residue_p, atom.padded_name()});
    }

    // Atom order from PDB: 0=N, 1=CA, 2=C, 3=O, 4=CB
    // Build 6 bonds mimicking what make_restraints produces
    // type 0 = BOND_RESTRAINT
    // N-CA, CA-C, CA-CB, CA-O, C-O, N-O
    restraints.restraints_vec.push_back({0, 0, 1});
    restraints.restraints_vec.push_back({0, 1, 2});
    restraints.restraints_vec.push_back({0, 1, 4});
    restraints.restraints_vec.push_back({0, 1, 3});
    restraints.restraints_vec.push_back({0, 2, 3});
    restraints.restraints_vec.push_back({0, 0, 3});

    // Case 1: delete a real bond (CA-C in ALA 21) - atom names are 4-char padded
    {
        size_t bonds_before = restraints.restraints_vec.size();
        EXPECT_EQ(bonds_before, 6);

        coot::chem_mod_bond mod_bond("delete", " CA ", " C  ", "", 0.0, 0.0);
        restraints.mod_bond_delete_gemmi(mod_bond, residue_p);

        size_t bonds_after = restraints.restraints_vec.size();
        EXPECT_EQ(bonds_after, 5);
    }

    // Case 2: try to delete a non-existent bond (atom names that don't exist)
    {
        size_t bonds_before2 = restraints.restraints_vec.size();
        EXPECT_EQ(bonds_before2, 5);

        coot::chem_mod_bond mod_bond("delete", "ZZZ ", "YYY ", "", 0.0, 0.0);
        restraints.mod_bond_delete_gemmi(mod_bond, residue_p);

        size_t bonds_after2 = restraints.restraints_vec.size();
        EXPECT_EQ(bonds_after2, 5);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}