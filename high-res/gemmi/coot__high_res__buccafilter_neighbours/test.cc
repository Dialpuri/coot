#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, buccafilter_neighbours) {
    // case: original oracle test
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        int original_count = 0;
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& res : chain.residues) {
                    original_count += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(original_count, 2107);
        
        coot::high_res::buccafilter_neighbours_gemmi(st);
        
        int final_count = 0;
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& res : chain.residues) {
                    final_count += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(final_count, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}