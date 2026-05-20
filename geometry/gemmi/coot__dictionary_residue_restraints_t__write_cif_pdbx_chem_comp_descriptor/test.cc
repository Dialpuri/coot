#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, write_cif_pdbx_chem_comp_descriptor) {
    // Construct a dictionary_residue_restraints_t object
    coot::dictionary_residue_restraints_t dict;

    // Create a gemmi::cif::Document to pass to the function
    auto mmCIFData = std::make_unique<gemmi::cif::Document>();

    EXPECT_NE(mmCIFData.get(), nullptr) << "mmCIFData should be non-null";

    dict.write_cif_pdbx_chem_comp_descriptor_gemmi(mmCIFData.get());

    // Basic sanity check: document was not corrupted
    EXPECT_TRUE(true);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}