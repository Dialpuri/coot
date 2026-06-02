#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, make_helix_pseudo_bond_restraints_from_res_vec_auto) {
    {
        // Case 1: Oracle's original case - residues 10-25 in chain A
        // This produces 0 restraints (helical_results.start == -1)
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_EQ(st.models.size(), 1u);

        gemmi::Model& model = st.models[0];

        // Collect consecutive residues from chain A (residues 10-25)
        std::vector<std::pair<bool, gemmi::Residue*>> residues;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                int seqnum = res.seqid.num.value;
                if (seqnum >= 10 && seqnum <= 25) {
                    residues.push_back(std::make_pair(false, &res));
                }
            }
        }

        EXPECT_EQ(residues.size(), 16u);

        int n_helical_restraints = coot::make_helix_pseudo_bond_restraints_from_res_vec_auto_gemmi(residues);

        // Original test: EXPECT_EQ(restraints_after, 0u)
        // and EXPECT_EQ((restraints_after - restraints_before), 0u)
        EXPECT_EQ(n_helical_restraints, 0);
    }

    {
        // Case 2: Try ALL protein residues in chain A to maximize chance of helix detection
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_EQ(st.models.size(), 1u);

        gemmi::Model& model = st.models[0];

        // Collect ALL residues from chain A
        std::vector<std::pair<bool, gemmi::Residue*>> residues;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                residues.push_back(std::make_pair(false, &res));
            }
        }

        EXPECT_FALSE(residues.empty());

        int n_helical_restraints = coot::make_helix_pseudo_bond_restraints_from_res_vec_auto_gemmi(residues);

        // Original test: EXPECT_GE(restraints_after, 0u) - just check it doesn't crash
        EXPECT_GE(n_helical_restraints, 0);
    }

    {
        // Case 3: Empty residue vector - should handle gracefully
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_EQ(st.models.size(), 1u);

        std::vector<std::pair<bool, gemmi::Residue*>> residues;

        EXPECT_NO_THROW(coot::make_helix_pseudo_bond_restraints_from_res_vec_auto_gemmi(residues));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}