#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

TEST(OracleTest, get_compound_lines) {
    // Case 1: PDB file (example.pdb) — oracle printed compound_lines_count: 0
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());
        std::vector<std::string> compound_lines = coot::get_compound_lines_gemmi(st);
        EXPECT_EQ(compound_lines.size(), 0u);
    }

    // Case 2: CIF file — example-ligand.cif
    {
        gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-ligand.cif");
        if (!st.models.empty()) {
            std::vector<std::string> compound_lines2 = coot::get_compound_lines_gemmi(st);
            EXPECT_EQ(compound_lines2.size(), 0u);
        }
    }

    // Case 3: protein-ligand CIF (oracle printed compound_lines_count: 0)
    {
        gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        EXPECT_FALSE(st.models.empty());
        std::vector<std::string> compound_lines3 = coot::get_compound_lines_gemmi(st);
        EXPECT_EQ(compound_lines3.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
