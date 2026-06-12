#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, fill_charged_atoms) {
    // Case 1: fill charges for first 3 atoms
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GT(st.models.size(), 0u);

        gemmi::Model &model = st.models[0];

        ASSERT_GT(model.chains.size(), 0u);
        gemmi::Chain &chain = model.chains[0];

        ASSERT_GT(chain.residues.size(), 0u);
        gemmi::Residue &residue = chain.residues[0];

        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(residue.seqid.num.value, 1);
        EXPECT_EQ(residue.name, "MET");

        EXPECT_EQ(residue.atoms.size(), 8u);

        // Build restraints with partial charges (matching real coot structure)
        coot::dictionary_residue_restraints_t rest;
        for (unsigned int i = 0; i < residue.atoms.size() && i < 3; i++) {
            coot::dict_atom da;
            da.atom_id = residue.atoms[i].name;
            da.atom_id_4c = residue.atoms[i].name;
            float charge_val = static_cast<float>(i) * 0.1f;
            da.partial_charge = std::make_pair(true, charge_val);
            rest.atom_info.push_back(da);
        }

        // Record charges BEFORE (gemmi charge is stored x100 as signed char)
        std::vector<int> before_charges;
        for (unsigned int i = 0; i < residue.atoms.size(); i++) {
            before_charges.push_back(residue.atoms[i].charge);
        }

        // Call fill_charged_atoms_gemmi
        coot::dipole::fill_charged_atoms_gemmi(&residue, rest);

        // Count changed atoms
        int changed_count = 0;
        for (unsigned int i = 0; i < residue.atoms.size(); i++) {
            int before = (i < before_charges.size()) ? before_charges[i] : 0;
            int after = residue.atoms[i].charge;
            if (before != after) {
                changed_count++;
            }
        }
        EXPECT_EQ(changed_count, 2);

        // Verify specific charges: CA=0.1, C=0.2 (gemmi charge is x100)
        for (unsigned int i = 0; i < residue.atoms.size(); i++) {
            gemmi::Atom &atom = residue.atoms[i];
            std::string atom_name(atom.name);
            if (atom_name == " CA ") {
                EXPECT_NEAR(static_cast<float>(atom.charge) / 100.0f, 0.1f, 1e-4);
            } else if (atom_name == " C ") {
                EXPECT_NEAR(static_cast<float>(atom.charge) / 100.0f, 0.2f, 1e-4);
            }
        }
    }

    // Case 2: empty restraints (no charges set)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GT(st.models.size(), 0u);

        gemmi::Model &model = st.models[0];

        ASSERT_GT(model.chains.size(), 0u);
        gemmi::Chain &chain = model.chains[0];

        gemmi::Residue *residue2 = nullptr;
        if (chain.residues.size() > 1) {
            residue2 = &chain.residues[1];
        } else {
            residue2 = &chain.residues[0];
        }
        ASSERT_NE(residue2, nullptr);

        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(residue2->seqid.num.value, 2);
        EXPECT_EQ(residue2->name, "GLU");

        coot::dictionary_residue_restraints_t rest2;
        // No partial_charges set

        // Record before
        std::vector<int> before2;
        for (unsigned int i = 0; i < residue2->atoms.size(); i++) {
            before2.push_back(residue2->atoms[i].charge);
        }

        coot::dipole::fill_charged_atoms_gemmi(residue2, rest2);

        // Check after
        bool any_charge_set = false;
        int after_count = 0;
        for (unsigned int i = 0; i < residue2->atoms.size(); i++) {
            int after = residue2->atoms[i].charge;
            int before = (i < (int)before2.size()) ? before2[i] : 0;
            if (before != after) {
                any_charge_set = true;
                after_count++;
            }
        }
        EXPECT_FALSE(any_charge_set);
        EXPECT_EQ(after_count, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
