#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, pdbcleanup_serial_residue_numbers) {
    // Load the test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Count residues before
    int nres_before = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            nres_before += chain.residues.size();
        }
    }

    EXPECT_EQ(nres_before, 268) << "Total residue count before cleanup";

    // Collect a sample of indices before (using group_idx as the storage field)
    std::vector<int> indices_before;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            int nres = chain.residues.size();
            for (int ires = 0; ires < std::min(5, nres); ++ires) {
                indices_before.push_back(chain.residues[ires].group_idx);
            }
        }
    }

    EXPECT_EQ(indices_before.size(), 6u) << "Sample indices count before cleanup";
    if (!indices_before.empty()) {
        EXPECT_EQ(indices_before[0], 0) << "First index before cleanup";
    }

    // Call the function
    coot::util::pdbcleanup_serial_residue_numbers_gemmi(st);

    // Check after: all indices should match their 0-based chain position
    bool all_indices_reset = true;
    int total_residues = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (int ires = 0; ires < static_cast<int>(chain.residues.size()); ++ires) {
                if (chain.residues[ires].group_idx != ires) {
                    all_indices_reset = false;
                }
                total_residues++;
            }
        }
    }

    EXPECT_TRUE(all_indices_reset) << "All residue indices should be reset to match 0-based chain position";
    EXPECT_EQ(total_residues, 268) << "Total residue count after cleanup";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}