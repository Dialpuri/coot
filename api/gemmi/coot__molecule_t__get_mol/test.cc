#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetMolAndNModels) {
   // Simulate molecules_container_t flow with direct molecule_t usage
   coot::molecule_t mol;
   mol.structure = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   gemmi::Structure* st = mol.get_mol_gemmi();
   EXPECT_NE(st, nullptr);
   
   if (st) {
      // Count models in gemmi - st.models.size() since gemmi::Structure has models vector
      size_t nModels = st->models.size();
      EXPECT_EQ(nModels, 1);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}