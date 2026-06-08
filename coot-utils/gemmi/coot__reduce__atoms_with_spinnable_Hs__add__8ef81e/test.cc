#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/api/gemmi/coot__molecule_t__cid_to_atom/function.hh"

TEST(OracleTest, atoms_with_spinnable_Hs_add) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    coot::reduce::atoms_with_spinnable_Hs spinables;

    // Case 1: valid atom with METHYL type and empty hydrogen vector
    {
        std::string cid = "//A/1/CA";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        std::vector<gemmi::Atom *> h_atoms;
        size_t before = spinables.typed_atoms[""].size();
        spinables.add_gemmi(at, coot::reduce::atom_with_attached_Hs::METHYL, h_atoms);
        size_t after = spinables.typed_atoms[""].size();

        EXPECT_NE(at, nullptr);
        EXPECT_EQ(h_atoms.size(), 0u);
        EXPECT_EQ(before, 0u);
        EXPECT_EQ(after, 1u);
    }

    // Case 2: valid atom with HYDROXYL type and one hydrogen atom in vector
    {
        std::string cid = "//A/2/CA";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        std::vector<gemmi::Atom *> h_atoms;
        // Reuse the atom itself as a stand-in hydrogen pointer
        if (at) h_atoms.push_back(at);
        size_t before = spinables.typed_atoms[""].size();
        spinables.add_gemmi(at, coot::reduce::atom_with_attached_Hs::HYDROXYL, h_atoms);
        size_t after = spinables.typed_atoms[""].size();

        EXPECT_NE(at, nullptr);
        EXPECT_EQ(h_atoms.size(), 1u);
        EXPECT_EQ(before, 1u);
        EXPECT_EQ(after, 2u);
    }

    // Case 3: atom with non-empty altLoc (tests alt_loc keying)
    {
        std::string cid = "//A/1/N";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        std::vector<gemmi::Atom *> h_atoms;
        std::string alt;
        if (at && at->altloc != ' ' && at->altloc != '\0') {
            alt = std::string(1, at->altloc);
        }
        size_t before = spinables.typed_atoms[alt].size();
        spinables.add_gemmi(at, coot::reduce::atom_with_attached_Hs::SULFHYDRYL, h_atoms);
        size_t after = spinables.typed_atoms[alt].size();

        EXPECT_NE(at, nullptr);
        EXPECT_EQ(h_atoms.size(), 0u);
        EXPECT_EQ(alt, "");
        EXPECT_EQ(before, 2u);
        EXPECT_EQ(after, 3u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
