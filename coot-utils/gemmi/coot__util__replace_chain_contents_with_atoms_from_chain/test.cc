#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <utility>
#include "function.hh"

TEST(OracleTest, replace_chain_contents_with_atoms_from_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    int imod = 0;  // gemmi is 0-indexed
    gemmi::Model& model = st.models[imod];
    ASSERT_FALSE(st.models.empty());

    // Pick chain A
    int ichain = 0;
    gemmi::Chain& chain = model.chains[ichain];
    ASSERT_FALSE(model.chains.empty());

    std::string chain_id(chain.name);
    EXPECT_EQ(chain_id, "A");

    int n_res_before = chain.residues.size();
    EXPECT_EQ(n_res_before, 267);

    // Case 1: Copy chain, trim it (remove last residue), then replace original contents
    {
        std::pair<gemmi::Chain*, gemmi::Structure> pair = coot::util::copy_chain_gemmi(chain);
        gemmi::Chain* copied_chain = pair.first;

        int n_res_copied = copied_chain->residues.size();
        EXPECT_EQ(n_res_copied, 267);

        // Remove the last residue from the copy
        if (n_res_copied > 1) {
            copied_chain->residues.erase(copied_chain->residues.end() - 1);
        }
        int n_res_after_trim = copied_chain->residues.size();
        EXPECT_EQ(n_res_after_trim, 266);

        bool do_finish = true;
        EXPECT_NO_THROW(
            coot::util::replace_chain_contents_with_atoms_from_chain_gemmi(&chain, &st, copied_chain, do_finish)
        );

        int n_res_after_replace = chain.residues.size();
        EXPECT_EQ(n_res_after_replace, 266);
        EXPECT_EQ(n_res_before - n_res_after_replace, 1);
    }

    // Case 2: do_finishstructedit = false — still copies, but doesn't call FinishStructEdit
    {
        std::pair<gemmi::Chain*, gemmi::Structure> pair = coot::util::copy_chain_gemmi(chain);
        gemmi::Chain* copied_chain = pair.first;

        // Delete all residues from the copy except the first
        int n_res_copied = copied_chain->residues.size();
        for (int i = n_res_copied - 1; i > 0; i--) {
            copied_chain->residues.erase(copied_chain->residues.begin() + i);
        }

        int n_res_in_copy = copied_chain->residues.size();
        int n_res_orig_before = chain.residues.size();
        EXPECT_EQ(n_res_in_copy, 1);
        EXPECT_EQ(n_res_orig_before, 266);

        bool do_finish = false;
        EXPECT_NO_THROW(
            coot::util::replace_chain_contents_with_atoms_from_chain_gemmi(&chain, &st, copied_chain, do_finish)
        );

        int n_res_orig_after = chain.residues.size();
        EXPECT_EQ(n_res_orig_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}