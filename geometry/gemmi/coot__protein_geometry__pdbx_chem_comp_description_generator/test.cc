#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, pdbx_chem_comp_description_generator) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Record initial state: count entries with non-empty description_generation
    int initial_count = 0;
    for (size_t i = 0; i < geom.size(); i++) {
        const auto &dg = geom[i].second.description_generation;
        if (!dg.program_name.empty() || !dg.program_version.empty() || !dg.descriptor.empty()) {
            initial_count++;
        }
    }
    EXPECT_EQ(initial_count, 34);
    EXPECT_EQ(geom.size(), 41u);

    // --- Case 1: valid comp_id (ALA) with description generator info ---
    {
        gemmi::cif::Loop loop;
        loop.tags.push_back("comp_id");
        loop.tags.push_back("program_name");
        loop.tags.push_back("program_version");
        loop.tags.push_back("descriptor");
        loop.values = {"ALA", "MyProgram", "1.0", "L-Alanine"};

        int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
        coot::pdbx_chem_comp_description_generator_gemmi(geom, loop, imol_enc);

        int alaidx = geom.get_monomer_restraints_index("ALA", imol_enc, true);
        bool found = (alaidx >= 0);
        EXPECT_TRUE(found);
        if (found) {
            const auto &dg = geom[alaidx].second.description_generation;
            EXPECT_EQ(dg.program_name, "MyProgram");
            EXPECT_EQ(dg.program_version, "1.0");
            EXPECT_EQ(dg.descriptor, "L-Alanine");
        }
    }

    // --- Case 2: invalid comp_id (ZZZ — not in dictionary) ---
    {
        gemmi::cif::Loop loop2;
        loop2.tags.push_back("comp_id");
        loop2.tags.push_back("program_name");
        loop2.tags.push_back("program_version");
        loop2.tags.push_back("descriptor");
        loop2.values = {"ZZZ", "FakeProgram", "0.1", "Nonexistent residue"};

        int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
        coot::pdbx_chem_comp_description_generator_gemmi(geom, loop2, imol_enc);

        int zzzidx = geom.get_monomer_restraints_index("ZZZ", imol_enc, true);
        bool found = (zzzidx >= 0);
        EXPECT_FALSE(found);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}