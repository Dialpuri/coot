#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue_types_without_dictionaries) {
    // Case 1: Standard PDB with init_standard() — all residues have dictionaries
    {
        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model &model = st.models[0];

        std::vector<std::string> result = coot::molecule_t::get_residue_types_without_dictionaries_gemmi(model, geom, 0);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: protein_geometry NOT initialized — no dictionaries loaded
    {
        // protein_geometry without init_standard() — no dictionary entries
        coot::protein_geometry geom;

        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model &model = st.models[0];

        std::vector<std::string> result = coot::molecule_t::get_residue_types_without_dictionaries_gemmi(model, geom, 0);

        // With no dictionary loaded, all protein-like residues should be reported
        EXPECT_GT(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}