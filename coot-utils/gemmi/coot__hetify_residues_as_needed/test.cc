#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: count atoms belonging to residues with het_flag == 'H'
int count_hetatms(const gemmi::Model& model) {
    int count = 0;
    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            if (res.het_flag == 'H') {
                count += static_cast<int>(res.atoms.size());
            }
        }
    }
    return count;
}

TEST(OracleTest, hetify_residues_as_needed) {
    // Case 1: standard protein with HOH in chain B
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int het_before = count_hetatms(st.models[0]);
        EXPECT_EQ(het_before, 1);

        int result = coot::hetify_residues_as_needed_gemmi(&st.models[0]);
        EXPECT_EQ(result, 0);

        int het_after = count_hetatms(st.models[0]);
        EXPECT_EQ(het_after, 1);
    }

    // Case 2: null pointer
    {
        int result = coot::hetify_residues_as_needed_gemmi(nullptr);
        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}