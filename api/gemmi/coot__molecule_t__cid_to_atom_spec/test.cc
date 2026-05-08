#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test for coot::molecule_t::cid_to_atom_spec_gemmi

TEST(OracleTest, CidToAtomSpec) {
   // Load structure using gemmi
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Try several CID formats to get valid inputs
   std::vector<std::string> cids = {
      "//A/10/CA",   // backbone atom of residue 10 in chain A
      "//A/1/CA",    // first residue
      "//A/2/O",     // oxygen of residue 2
      "//A/3/CB",    // beta carbon of residue 3
   };

   for (const auto &cid : cids) {
      std::pair<bool, coot::atom_spec_t> result = coot::cid_to_atom_spec_gemmi(st, cid);
      EXPECT_TRUE(result.first);
      if (result.first) {
         const coot::atom_spec_t &spec = result.second;
         EXPECT_EQ(spec.chain_id, "A");
         // Verify residue number matches the CID
         int expected_res_no = (cid == "//A/10/CA") ? 10 : 
                              (cid == "//A/1/CA") ? 1 :
                              (cid == "//A/2/O") ? 2 : 3;
         EXPECT_EQ(spec.res_no, expected_res_no);
         // Verify atom name - left-padded to 4 characters
         std::string expected_atom = (cid == "//A/10/CA") ? " CA " :
                                    (cid == "//A/1/CA") ? " CA " :
                                    (cid == "//A/2/O") ? " O  " : " CB ";
         EXPECT_EQ(spec.atom_name, expected_atom);
         // ins_code should be empty
         EXPECT_EQ(spec.ins_code, "");
         // alt_conf should be empty
         EXPECT_EQ(spec.alt_conf, "");
      }
   }

   // Edge case: invalid CID
   std::pair<bool, coot::atom_spec_t> result = coot::cid_to_atom_spec_gemmi(st, "//Z/999/CA");
   EXPECT_FALSE(result.first);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}