#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ChangeChainIdInsertOrAdd) {
  // Load test PDB file
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  
  // Get the first model
  ASSERT_GE(st.models.size(), 1);
  gemmi::Model& model = st.models[0];
  
  // Get the first chain (chain A)
  ASSERT_GE(model.chains.size(), 1);
  gemmi::Chain& chain_a = model.chains[0];
  
  printf("Model has %zu chains\n", model.chains.size());
  printf("Chain A has %zu residues\n", chain_a.residues.size());
  
  // Test 2: Insert in the middle (InsertResidue path)
  {
    printf("Creating chain_c...\n");
    gemmi::Chain chain_c("C");
    
    // Manually add a residue to chain_c before pushing
    gemmi::Residue res;
    res.seqid = chain_a.residues[19].seqid;
    res.name = chain_a.residues[19].name;
    res.atoms.push_back(chain_a.residues[19].atoms[0]);
    chain_c.residues.push_back(res);
    
    printf("chain_c before push_back: %zu residues\n", chain_c.residues.size());
    
    model.chains.push_back(chain_c);
    printf("After push_back, model.chains.size() = %zu\n", model.chains.size());
    
    gemmi::Chain& chain_c_ref = model.chains.back();
    printf("chain_c_ref after push_back: %zu residues\n", chain_c_ref.residues.size());
    
    EXPECT_EQ(chain_c_ref.residues.size(), 1);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
