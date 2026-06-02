#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, stack_and_pair) {
    // Create and load PDB
    // In the MMDB original, mol->ReadCoorFile() returns false on success (counterintuitive),
    // so EXPECT_FALSE(loaded) means "the file loaded successfully".
    // In gemmi, read_pdb_file throws on error, so if we get here, loading succeeded.
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    EXPECT_EQ(st.models.size(), 1);
    
    // Original test exits early when !loaded, so stack_and_pair cases are never reached.
    // In gemmi, if we got here, loading succeeded, so we can now call the function.
    gemmi::Model& model = st.models[0];
    auto normal_map = coot::stack_and_pair::stack_and_pair_gemmi(model);
    // The original test doesn't assert on the map contents, so we just verify it was created.
    // (The map may contain 0 or more entries depending on the residues in the file.)
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}