#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(SimpleGeometricDistortions, FromMol) {
   // Load the test structure
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Get the first model
   gemmi::Model &model = st.models[0];
   
   // Find residue ILE 10 (sequence number 10)
   gemmi::CRA cra_found{nullptr, nullptr, nullptr};
   for (auto& chain : model.chains) {
      for (auto& res : chain.residues) {
         if (res.seqid.num.value == 10) {
            cra_found = {&chain, &res, nullptr};
            break;
         }
      }
      if (cra_found.residue) break;
   }
   
   ASSERT_NE(cra_found.residue, nullptr) << "Residue ILE 10 not found";
   
   // Test case 1: with_nbcs = false
   {
      bool with_nbcs = false;
      std::pair<int, double> result = coot::simple_geometric_distortions_from_mol_gemmi(cra_found, with_nbcs);
      EXPECT_EQ(result.first, 1);
      EXPECT_NEAR(result.second, 47.8626, 0.0001);
   }
   
   // Test case 2: with_nbcs = true
   {
      bool with_nbcs = true;
      std::pair<int, double> result = coot::simple_geometric_distortions_from_mol_gemmi(cra_found, with_nbcs);
      EXPECT_EQ(result.first, 1);
      EXPECT_NEAR(result.second, 39.3672, 0.0001);
   }
}

TEST(SimpleGeometricDistortions, InvalidResidue) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   gemmi::Model &model = st.models[0];
   
   // Find a non-existent residue (seqid 999)
   gemmi::CRA cra{nullptr, nullptr, nullptr};
   for (auto& chain : model.chains) {
      for (auto& res : chain.residues) {
         if (res.seqid.num.value == 999) {
            cra = {&chain, &res, nullptr};
            break;
         }
      }
      if (cra.residue) break;
   }
   
   // If no such residue exists, pass a CRA with a null residue
   std::pair<int, double> result = coot::simple_geometric_distortions_from_mol_gemmi(cra, false);
   EXPECT_EQ(result.first, 0);
   EXPECT_FLOAT_EQ(result.second, 0.0);
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}