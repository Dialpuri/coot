#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

TEST(OracleTest, new_positions_for_residue_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get model 0 (0-indexed in gemmi; MMDB was 1-indexed)
    gemmi::Model& model = st.models[0];

    // Get a valid residue pointer for //A/1 (MET)
    gemmi::Residue* residue_p = coot::cid_to_residue_gemmi("//A/1", model);
    EXPECT_NE(residue_p, nullptr);

    // --- Case 1: Valid residue, matching atom name (CA) ---
    {
        std::vector<coot::api::moved_atom_t> moved_atoms;
        moved_atoms.push_back(coot::api::moved_atom_t(" CA ", "", 100.0f, 200.0f, 300.0f));

        int n_moved = coot::new_positions_for_residue_atoms_gemmi(residue_p, moved_atoms, false);

        EXPECT_EQ(n_moved, 1);

        // Verify position actually changed
        for (auto& atom : residue_p->atoms) {
            if (atom.name == "CA") {
                EXPECT_NEAR(atom.pos.x, 100.0, 1e-4);
                EXPECT_NEAR(atom.pos.y, 200.0, 1e-4);
                EXPECT_NEAR(atom.pos.z, 300.0, 1e-4);
                break;
            }
        }
    }

    // --- Case 2: Valid residue, non-matching atom name ---
    {
        std::vector<coot::api::moved_atom_t> moved_atoms;
        moved_atoms.push_back(coot::api::moved_atom_t("ZZZ", "", 0.0f, 0.0f, 0.0f));

        int n_moved = coot::new_positions_for_residue_atoms_gemmi(residue_p, moved_atoms, false);

        EXPECT_EQ(n_moved, 0);
    }

    // --- Case 3: Null residue pointer ---
    {
        std::vector<coot::api::moved_atom_t> moved_atoms;
        moved_atoms.push_back(coot::api::moved_atom_t("CA", "", 0.0f, 0.0f, 0.0f));

        int n_moved = coot::new_positions_for_residue_atoms_gemmi(nullptr, moved_atoms, false);

        EXPECT_EQ(n_moved, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}