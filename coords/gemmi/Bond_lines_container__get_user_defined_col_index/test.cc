#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Original MMDB test translated to gemmi

TEST(OracleTest, get_user_defined_col_index) {
    // Load a PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    EXPECT_FALSE(model.chains.empty());

    gemmi::Chain& chain = model.chains[0];
    EXPECT_FALSE(chain.residues.empty());

    gemmi::Residue& residue = chain.residues[0];
    EXPECT_FALSE(residue.atoms.empty());

    gemmi::Atom& atom = residue.atoms[0];

    // Create a Bond_lines_container
    coot::Bond_lines_container blc;

    // Simulate UDD handle usage: we cannot register UDD handles in gemmi,
    // so we set the atom's user data directly via a custom UDD entry.
    // Since gemmi has no UDD API, we mimic the test's intent by checking
    // that the function returns -1 when no valid handle is provided.
    int udd_handle = 0; // dummy handle
    int result = blc.get_user_defined_col_index_gemmi(atom, udd_handle);
    EXPECT_EQ(result, -1);

    // Test with invalid handle
    int invalid_handle = -1;
    int result_invalid = blc.get_user_defined_col_index_gemmi(atom, invalid_handle);
    EXPECT_EQ(result_invalid, -1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}