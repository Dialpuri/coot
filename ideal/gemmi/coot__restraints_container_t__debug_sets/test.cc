#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_t_debug_sets) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());

    // Case 1: Single residue (moving) - first residue in chain
    {
        gemmi::Residue* residue_p = &(chain.residues[0]);
        ASSERT_NE(residue_p, nullptr);

        coot::restraints_container_t restraints;
        restraints.mol = &model;
        restraints.residues_vec.push_back(std::make_pair(false, residue_p));

        // debug_sets_gemmi prints diagnostic info to cout - should not throw
        EXPECT_NO_THROW(restraints.debug_sets_gemmi());
    }

    // Case 2: Two residues (moving + fixed) - first two residues in chain
    {
        ASSERT_GE(chain.residues.size(), 2u);
        gemmi::Residue* res1 = &(chain.residues[0]); // first residue
        gemmi::Residue* res2 = &(chain.residues[1]); // second residue
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);

        coot::restraints_container_t restraints2;
        restraints2.mol = &model;
        restraints2.residues_vec.push_back(std::make_pair(false, res1)); // moving
        restraints2.residues_vec.push_back(std::make_pair(true, res2));  // fixed

        // Add fixed neighbours to test that code path too
        restraints2.fixed_neighbours_set[res1].insert(res2);

        // debug_sets_gemmi prints diagnostic info to cout - should not throw
        EXPECT_NO_THROW(restraints2.debug_sets_gemmi());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}