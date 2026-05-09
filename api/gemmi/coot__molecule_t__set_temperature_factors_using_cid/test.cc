#define _NO_PYTHON
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, set_temperature_factors_using_cid) {
   // Test 1: Valid CID //A/1
   {
      auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      // Get initial tempFactor for //A/1
      float initial_temp = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 1) {
              if (!res.atoms.empty())
                initial_temp = res.atoms[0].b_iso;
              break;
            }
          }
      
      // Call the function
      coot::molecule_t::set_temperature_factors_using_cid_gemmi(st, "//A/1", 50.0f);
      
      // Verify change
      float final_temp = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 1) {
              if (!res.atoms.empty())
                final_temp = res.atoms[0].b_iso;
              break;
            }
          }
      
      EXPECT_FLOAT_EQ(final_temp, 50.0f);
   }
   
   // Test 2: Valid CID //A/10
   {
      auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      // Get initial tempFactor for //A/10
      float initial_temp_10 = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 10) {
              if (!res.atoms.empty())
                initial_temp_10 = res.atoms[0].b_iso;
              break;
            }
          }
      
      // Call the function
      coot::molecule_t::set_temperature_factors_using_cid_gemmi(st, "//A/10", 75.0f);
      
      // Verify change
      float final_temp_10 = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 10) {
              if (!res.atoms.empty())
                final_temp_10 = res.atoms[0].b_iso;
              break;
            }
          }
      
      EXPECT_FLOAT_EQ(final_temp_10, 75.0f);
   }
   
   // Test 3: Invalid CID //X/999 (should not crash, should not modify anything)
   {
      auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      // Get initial tempFactor for //A/1 (before invalid CID)
      float temp_before_invalid = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 1) {
              if (!res.atoms.empty())
                temp_before_invalid = res.atoms[0].b_iso;
              break;
            }
          }
      
      // Call with invalid CID (should silently do nothing)
      coot::molecule_t::set_temperature_factors_using_cid_gemmi(st, "//X/999", 99.0f);
      
      // Verify //A/1 was not modified
      float temp_after_invalid = -999.0f;
      for (auto& model : st.models)
        for (auto& chain : model.chains)
          for (auto& res : chain.residues) {
            if (chain.name == "A" && res.seqid.num.value == 1) {
              if (!res.atoms.empty())
                temp_after_invalid = res.atoms[0].b_iso;
              break;
            }
          }
      
      EXPECT_FLOAT_EQ(temp_after_invalid, temp_before_invalid);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}