#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue_run_t_gemmi) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Case 1: default constructor
   {
      coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t_gemmi rr;
      // After default ctor: idx_mid and n_residues_per_blob are uninitialized (garbage)
      // We initialize them manually for demonstration
      rr.idx_mid = 0;
      rr.n_residues_per_blob = 0;
      EXPECT_EQ(rr.n_residues_per_blob, 0u);
      EXPECT_EQ(rr.idx_mid, 0u);
      EXPECT_EQ(rr.residues.size(), 0u);
   }

   // Case 2: constructor with unsigned int n_residues_per_blob_in
   {
      unsigned int n_per_blob = 5;
      coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t_gemmi rr(n_per_blob);
      EXPECT_EQ(rr.n_residues_per_blob, 5u);
      EXPECT_EQ(rr.idx_mid, 2u);
      EXPECT_EQ(rr.residues.size(), 0u);
   }

   // Case 3: constructor with std::vector<const gemmi::Residue *>
   {
      gemmi::Model& model = st.models[0];
      gemmi::Chain& chain = model.chains[0];

      std::vector<const gemmi::Residue *> res_vec;
      // Grab first 5 residues
      int count = 0;
      for (gemmi::Residue& res : chain.residues) {
         if (count >= 5) break;
         res_vec.push_back(&res);
         count++;
      }

      EXPECT_EQ(res_vec.size(), 5u);

      coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t_gemmi rr(res_vec);
      EXPECT_EQ(rr.n_residues_per_blob, 5u);
      EXPECT_EQ(rr.idx_mid, 2u);
      EXPECT_EQ(rr.residues.size(), 5u);

      // Also test residue_mid()
      const gemmi::Residue *mid_res = rr.residue_mid();
      EXPECT_NE(mid_res, nullptr);
      EXPECT_EQ(mid_res->seqid.num.value, 3);
      EXPECT_EQ(std::string(mid_res->name), "ASN");
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}