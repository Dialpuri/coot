#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dots_representation_info_t) {
    // Load the PDB - try multiple file paths
    gemmi::Structure st;
    bool loaded = false;

    // Try example.pdb first
    try {
        st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        loaded = true;
    } catch (...) {}

    if (!loaded) {
        try {
            st = gemmi::read_pdb_file("example.pdb");
            loaded = true;
        } catch (...) {}
    }

    if (!loaded) {
        try {
            st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
            loaded = true;
        } catch (...) {}
    }

    std::cout << "loaded: " << loaded << " models: " << st.models.size() << std::endl;
    ASSERT_TRUE(loaded);

    gemmi::Structure const* mol = &st;

    // Case 1: Construct with mol and mol_exclude = nullptr
    {
        gemmi::Structure const* mol_exclude = nullptr;
        pli::dots_representation_info_t dots = pli::dots_representation_info_t_gemmi(mol, mol_exclude);
        bool is_open_val = dots.is_open_p();
        std::cout << "Case 1 is_open: " << (is_open_val ? "true" : "false") << std::endl;
        EXPECT_TRUE(is_open_val);
    }

    // Case 2: Construct with mol and mol_exclude = same mol (self-exclusion)
    {
        pli::dots_representation_info_t dots = pli::dots_representation_info_t_gemmi(mol, mol);
        bool is_open_val = dots.is_open_p();
        std::cout << "Case 2 is_open: " << (is_open_val ? "true" : "false") << std::endl;
        EXPECT_TRUE(is_open_val);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}