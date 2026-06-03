#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_comp_synonym) {
    coot::protein_geometry::residue_name_synonym ala("ALA", "ALA_OLD", "default");
    coot::protein_geometry::residue_name_synonym gly("GLY", "GLY_ALT", "default");
    coot::protein_geometry::residue_name_synonym ser("SER", "SER_X", "");
    
    std::vector<coot::protein_geometry::residue_name_synonym> residue_name_synonyms;
    residue_name_synonyms.push_back(ala);
    residue_name_synonyms.push_back(gly);
    residue_name_synonyms.push_back(ser);
    
    size_t initial_size = residue_name_synonyms.size();
    EXPECT_EQ(initial_size, 3U);

    // Case 1: main oracle case — adds 3 synonym entries (success path)
    {
        gemmi::cif::Loop loop;
        loop.tags.push_back("comp_id");
        loop.tags.push_back("comp_alternative_id");
        loop.tags.push_back("mod_id");

        loop.values.push_back("ALA");
        loop.values.push_back("ALA_OLD");
        loop.values.push_back("default");

        loop.values.push_back("GLY");
        loop.values.push_back("GLY_ALT");
        loop.values.push_back("default");

        loop.values.push_back("SER");
        loop.values.push_back("SER_X");
        loop.values.push_back("");

        EXPECT_EQ(loop.length(), 3U);
        EXPECT_EQ(loop.val(0, 0), "ALA");
        EXPECT_EQ(loop.val(0, 1), "ALA_OLD");
        EXPECT_EQ(loop.val(0, 2), "default");

        EXPECT_NO_THROW(coot::protein_geometry::add_chem_comp_synonym_gemmi(loop, residue_name_synonyms));

        EXPECT_EQ(residue_name_synonyms.size(), initial_size + 3);

        size_t n = initial_size;
        EXPECT_EQ(residue_name_synonyms[n].comp_id, "ALA");
        EXPECT_EQ(residue_name_synonyms[n].comp_alternative_id, "ALA_OLD");
        EXPECT_EQ(residue_name_synonyms[n].mod_id, "default");
        n++;

        EXPECT_EQ(residue_name_synonyms[n].comp_id, "GLY");
        EXPECT_EQ(residue_name_synonyms[n].comp_alternative_id, "GLY_ALT");
        EXPECT_EQ(residue_name_synonyms[n].mod_id, "default");
        n++;

        EXPECT_EQ(residue_name_synonyms[n].comp_id, "SER");
        EXPECT_EQ(residue_name_synonyms[n].comp_alternative_id, "SER_X");
        EXPECT_EQ(residue_name_synonyms[n].mod_id, "");
    }

    // Case 2: coverage compensation — empty loop (no rows, exercises early return branch)
    {
        size_t size_before = residue_name_synonyms.size();

        gemmi::cif::Loop loop;
        loop.tags.push_back("comp_id");
        loop.tags.push_back("comp_alternative_id");
        loop.tags.push_back("mod_id");
        // No rows added

        EXPECT_EQ(loop.length(), 0U);

        EXPECT_NO_THROW(coot::protein_geometry::add_chem_comp_synonym_gemmi(loop, residue_name_synonyms));

        EXPECT_EQ(residue_name_synonyms.size(), size_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}