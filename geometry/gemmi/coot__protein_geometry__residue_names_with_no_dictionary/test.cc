#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, residue_names_with_no_dictionary) {
    // Initialize protein geometry with standard dictionaries
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Standard protein with mostly standard residues + HOH water
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int imol_no = 0;

        std::vector<std::string> result = coot::residue_names_with_no_dictionary_gemmi(geom, &st, imol_no);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: Null pointer — guarded path, should return empty
    {
        std::vector<std::string> result = coot::residue_names_with_no_dictionary_gemmi(geom, nullptr, 0);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: Non-standard residue name (complementary case)
    // Create a structure with a residue that has no dictionary entry
    {
        gemmi::Structure st;
        gemmi::Model model;
        gemmi::Chain chain;
        chain.name = "A";

        gemmi::Residue residue;
        residue.name = "XFOO";
        residue.seqid = gemmi::SeqId(1, ' ');

        gemmi::Atom atom;
        atom.name = "CA";
        atom.pos = gemmi::Position(0, 0, 0);
        atom.element = gemmi::Element("C");
        atom.altloc = '\0';
        atom.occ = 1.0;
        atom.b_iso = 0.0;
        residue.atoms.push_back(atom);

        chain.residues.push_back(residue);
        model.chains.push_back(chain);
        st.models.push_back(model);

        int imol_no = 0;

        std::vector<std::string> result = coot::residue_names_with_no_dictionary_gemmi(geom, &st, imol_no);

        EXPECT_EQ(result.size(), 1u);
        EXPECT_EQ(result[0], "XFOO");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}