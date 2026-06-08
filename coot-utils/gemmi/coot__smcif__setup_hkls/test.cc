#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Case 1: protein-ligand CIF (likely no refln index loops)
TEST(OracleTestGemmi, setup_hkls_protein_ligand) {
    std::string cif = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif";
    
    auto hkls = coot::smcif::setup_hkls_gemmi(cif);
    
    EXPECT_EQ(hkls.size(), 0);
}

// Case 2: nonexistent file — tests the error guard path
TEST(OracleTestGemmi, setup_hkls_nonexistent) {
    std::string cif = "/nonexistent.cif";
    
    auto hkls = coot::smcif::setup_hkls_gemmi(cif);
    
    EXPECT_EQ(hkls.size(), 0);
}

// Case 3: CIF with actual refln index loop — exercise the success path
TEST(OracleTestGemmi, setup_hkls_with_refln) {
    const char* tmp_cif = "/tmp/test_refln_hkls.cif";
    FILE* f = fopen(tmp_cif, "w");
    if (!f) {
        GTEST_SKIP() << "Cannot write temp CIF file";
        return;
    }
    fprintf(f, "data_test\n");
    fprintf(f, "loop_\n");
    fprintf(f, "_refln_index_h\n");
    fprintf(f, "_refln_index_k\n");
    fprintf(f, "_refln_index_l\n");
    fprintf(f, "1 2 3\n");
    fprintf(f, "4 5 6\n");
    fprintf(f, "7 8 9\n");
    fclose(f);
    
    auto hkls = coot::smcif::setup_hkls_gemmi(tmp_cif);
    
    EXPECT_EQ(hkls.size(), 3);
    
    // Check individual values using accessor methods
    EXPECT_EQ(hkls[0].h(), 1);
    EXPECT_EQ(hkls[0].k(), 2);
    EXPECT_EQ(hkls[0].l(), 3);
    EXPECT_EQ(hkls[1].h(), 4);
    EXPECT_EQ(hkls[1].k(), 5);
    EXPECT_EQ(hkls[1].l(), 6);
    EXPECT_EQ(hkls[2].h(), 7);
    EXPECT_EQ(hkls[2].k(), 8);
    EXPECT_EQ(hkls[2].l(), 9);
    
    unlink(tmp_cif);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}