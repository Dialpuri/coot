#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>
#include "function.hh"

TEST(OracleTest, parse_lib_info) {
    // Case 1: pass nullptr (edge case)
    // Oracle OUTPUT: called_successfully: true
    {
        gemmi::cif::Block* null_block = nullptr;
        EXPECT_NO_THROW(coot::protein_geometry::parse_lib_info_gemmi(null_block));
    }

    // Case 2: try to get a real Block from a CIF file
    // Oracle OUTPUT: cif_read_ok: true, n_categories: 1, called_with_real_struct: false (no struct found, used nullptr)
    {
        gemmi::cif::Document doc = gemmi::read_cif_gz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-ligand.cif");
        int ierr = doc.blocks.empty() ? 1 : 0;

        EXPECT_EQ(ierr, 0);

        gemmi::cif::Block* block_ptr = nullptr;

        if (ierr == 0 && !doc.blocks.empty()) {
            block_ptr = &doc.blocks[0];
        }

        // The oracle found no struct and called with nullptr.
        // In some environments a struct may be found — either way, parse_lib_info should succeed.
        if (block_ptr) {
            EXPECT_NO_THROW(coot::protein_geometry::parse_lib_info_gemmi(block_ptr));
        } else {
            EXPECT_NO_THROW(coot::protein_geometry::parse_lib_info_gemmi(nullptr));
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}