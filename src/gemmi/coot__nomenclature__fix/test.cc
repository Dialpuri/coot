#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "function.hh"

TEST(OracleTest, nomenclature_fix_and_list) {
    // Load PDB via gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Build protein geometry with standard dictionary
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: call fix on the full model
    {
        std::vector<gemmi::CRA> fixed = coot::nomenclature::fix_gemmi(st.models[0], &geom);

        EXPECT_EQ(fixed.size(), 8u);

        // Verify each fixed residue matches oracle output
        EXPECT_NE(fixed[0].residue, nullptr);
        EXPECT_EQ(std::string(fixed[0].chain->name), "A");
        EXPECT_EQ(fixed[0].residue->seqid.num.value, 4);
        EXPECT_EQ(std::string(fixed[0].residue->name), "PHE");

        EXPECT_NE(fixed[1].residue, nullptr);
        EXPECT_EQ(std::string(fixed[1].chain->name), "A");
        EXPECT_EQ(fixed[1].residue->seqid.num.value, 8);
        EXPECT_EQ(std::string(fixed[1].residue->name), "GLU");

        EXPECT_NE(fixed[2].residue, nullptr);
        EXPECT_EQ(std::string(fixed[2].chain->name), "A");
        EXPECT_EQ(fixed[2].residue->seqid.num.value, 19);
        EXPECT_EQ(std::string(fixed[2].residue->name), "TYR");

        EXPECT_NE(fixed[3].residue, nullptr);
        EXPECT_EQ(std::string(fixed[3].chain->name), "A");
        EXPECT_EQ(fixed[3].residue->seqid.num.value, 28);
        EXPECT_EQ(std::string(fixed[3].residue->name), "GLU");

        EXPECT_NE(fixed[4].residue, nullptr);
        EXPECT_EQ(std::string(fixed[4].chain->name), "A");
        EXPECT_EQ(fixed[4].residue->seqid.num.value, 82);
        EXPECT_EQ(std::string(fixed[4].residue->name), "PHE");

        EXPECT_NE(fixed[5].residue, nullptr);
        EXPECT_EQ(std::string(fixed[5].chain->name), "A");
        EXPECT_EQ(fixed[5].residue->seqid.num.value, 208);
        EXPECT_EQ(std::string(fixed[5].residue->name), "GLU");

        EXPECT_NE(fixed[6].residue, nullptr);
        EXPECT_EQ(std::string(fixed[6].chain->name), "A");
        EXPECT_EQ(fixed[6].residue->seqid.num.value, 236);
        EXPECT_EQ(std::string(fixed[6].residue->name), "TYR");

        EXPECT_NE(fixed[7].residue, nullptr);
        EXPECT_EQ(std::string(fixed[7].chain->name), "A");
        EXPECT_EQ(fixed[7].residue->seqid.num.value, 247);
        EXPECT_EQ(std::string(fixed[7].residue->name), "ASP");
    }

    // Case 2: re-call fix on already-fixed model — should return 0
    {
        std::vector<gemmi::CRA> listed = coot::nomenclature::fix_gemmi(st.models[0], &geom);
        EXPECT_EQ(listed.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}