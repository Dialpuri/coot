#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomSpecMatchesSpec) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model& model = st.models[0];
  
  gemmi::Chain* chain_a = nullptr;
  for (auto& chain : model.chains) {
    if (chain.name == "A") {
      chain_a = &chain;
      break;
    }
  }
  ASSERT_TRUE(chain_a != nullptr);
  
  // Get residue 10 (ILE)
  gemmi::Residue* residue_10 = nullptr;
  for (auto& res : chain_a->residues) {
    if (res.seqid.num.value == 10) {
      residue_10 = &res;
      break;
    }
  }
  ASSERT_TRUE(residue_10 != nullptr);
  
  // Get first atom (N)
  ASSERT_TRUE(!residue_10->atoms.empty());
  gemmi::Atom* atom_n = &residue_10->atoms[0];
  gemmi::CRA cra{chain_a, residue_10, atom_n};
  
  // Test Case 1: Match first atom (N) completely
  {
    coot::atom_spec_t spec;
    spec.chain_id = cra.chain->name;
    spec.res_no = cra.residue->seqid.num.value;
    spec.ins_code = std::string(1, cra.residue->seqid.icode);
    spec.atom_name = cra.atom->name;
    spec.alt_conf = std::string(1, cra.atom->altloc);
    
    bool result = spec.matches_spec_gemmi(cra);
    EXPECT_TRUE(result);
  }
  
  // Test Case 2: Mismatch atom name (N -> CA)
  {
    coot::atom_spec_t spec;
    spec.chain_id = cra.chain->name;
    spec.res_no = cra.residue->seqid.num.value;
    spec.ins_code = std::string(1, cra.residue->seqid.icode);
    spec.atom_name = "CA"; // different atom name
    spec.alt_conf = std::string(1, cra.atom->altloc);
    
    bool result = spec.matches_spec_gemmi(cra);
    EXPECT_FALSE(result);
  }
  
  // Test Case 3: Mismatch chain ID (A -> B)
  {
    coot::atom_spec_t spec;
    spec.chain_id = "B"; // different chain
    spec.res_no = cra.residue->seqid.num.value;
    spec.ins_code = std::string(1, cra.residue->seqid.icode);
    spec.atom_name = cra.atom->name;
    spec.alt_conf = std::string(1, cra.atom->altloc);
    
    bool result = spec.matches_spec_gemmi(cra);
    EXPECT_FALSE(result);
  }
  
  // Test Case 4: Match valid atom with residue 21 (LYS)
  gemmi::Residue* residue_21 = nullptr;
  for (auto& res : chain_a->residues) {
    if (res.seqid.num.value == 21) {
      residue_21 = &res;
      break;
    }
  }
  if (residue_21 != nullptr && !residue_21->atoms.empty()) {
    gemmi::Atom* atom = &residue_21->atoms[0];
    gemmi::CRA cra21{chain_a, residue_21, atom};
    
    coot::atom_spec_t spec;
    spec.chain_id = "A";
    spec.res_no = cra21.residue->seqid.num.value;
    spec.ins_code = std::string(1, cra21.residue->seqid.icode);
    spec.atom_name = atom->name;
    spec.alt_conf = std::string(1, atom->altloc);
    
    bool result = spec.matches_spec_gemmi(cra21);
    EXPECT_TRUE(result);
  }
  
  // Test Case 5: Non-existent residue (res_no 999)
  {
    coot::atom_spec_t spec;
    spec.chain_id = "A";
    spec.res_no = 999;
    spec.ins_code = " ";
    spec.atom_name = "CA";
    spec.alt_conf = " ";
    
    bool result = spec.matches_spec_gemmi(cra);
    EXPECT_FALSE(result);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}