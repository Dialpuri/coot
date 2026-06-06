#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

// Helper: create a Structure containing only the specified chains from a model
static gemmi::Structure filter_chains(const gemmi::Model& model, std::vector<std::string> chain_names) {
   gemmi::Structure st;
   gemmi::Model new_model(1);
   st.models.push_back(std::move(new_model));
   gemmi::Model& dest_model = st.models[0];
   for (const auto& chain : model.chains) {
      bool keep = false;
      for (const auto& cn : chain_names) {
         if (chain.name == cn) { keep = true; break; }
      }
      if (keep) {
         gemmi::Chain new_chain;
         new_chain.name = chain.name;
         for (const auto& res : chain.residues) {
            new_chain.residues.push_back(res);
         }
         dest_model.chains.push_back(std::move(new_chain));
      }
   }
   return st;
}

TEST(OracleTest, insert_coords) {
    // Case 1: Original oracle case — insert all atoms from example.pdb
    {
        gemmi::Structure target = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = target.models[0];
        int n_chains_before = static_cast<int>(model.chains.size());
        int n_res_before = 0;
        if (!model.chains.empty()) {
            n_res_before = static_cast<int>(model.chains[0].residues.size());
        }

        EXPECT_EQ(n_chains_before, 2);
        EXPECT_EQ(n_res_before, 267);

        // Create a second structure (the fragment to insert)
        gemmi::Structure frag_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& frag_model = frag_st.models[0];

        // Count selected atoms (all atoms)
        int n_sel = static_cast<int>(gemmi::count_atom_sites(frag_model));
        EXPECT_EQ(n_sel, 2107);

        // Call insert_coords
        EXPECT_NO_THROW(insert_coords_gemmi(target, frag_model));

        // Check results
        gemmi::Model& model_after = target.models[0];
        int n_chains_after = static_cast<int>(model_after.chains.size());
        int n_res_after = 0;
        if (!model_after.chains.empty()) {
            n_res_after = static_cast<int>(model_after.chains[0].residues.size());
        }

        EXPECT_EQ(n_chains_after, 2);
        EXPECT_EQ(n_res_after, 534);

        // Complementary check: residues added = 534 - 267 = 267
        EXPECT_EQ(n_res_after - n_res_before, 267);
    }

    // Case 2: Select only atoms from chain A
    // Insert only a subset of atoms — exercises the loop with a smaller selection
    {
        gemmi::Structure target = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = target.models[0];
        int n_chains_before = static_cast<int>(model.chains.size());
        int n_res_before = 0;
        if (!model.chains.empty()) {
            n_res_before = static_cast<int>(model.chains[0].residues.size());
        }

        // Create a second structure (the fragment to insert)
        gemmi::Structure frag_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Select only atoms from chain A by creating a filtered structure
        gemmi::Structure frag_filtered = filter_chains(frag_st.models[0], {"A"});
        gemmi::Model& frag_model = frag_filtered.models[0];

        int n_sel = static_cast<int>(gemmi::count_atom_sites(frag_model));
        EXPECT_GT(n_sel, 0);

        // Call insert_coords with only chain A atoms
        EXPECT_NO_THROW(insert_coords_gemmi(target, frag_model));

        // Check results — chain count should still be 2
        gemmi::Model& model_after = target.models[0];
        int n_chains_after = static_cast<int>(model_after.chains.size());
        int n_res_after = 0;
        if (!model_after.chains.empty()) {
            n_res_after = static_cast<int>(model_after.chains[0].residues.size());
        }

        EXPECT_EQ(n_chains_after, 2);
        // Residues should increase (at least by 1 since we inserted atoms forming new residues)
        EXPECT_GT(n_res_after, n_res_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}