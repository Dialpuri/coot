#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string trim_atom_name(const std::string &s) {
    std::string r = s;
    // trim leading spaces
    auto start = r.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    r = r.substr(start);
    // trim trailing spaces
    while (!r.empty() && r.back() == ' ') r.pop_back();
    return r;
}

TEST(OracleTest, cid_to_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Test case 1: Valid atom CID (chain A, residue 10, atom N)
    {
        std::string cid = "//A/10/N";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(trim_atom_name(atom->name), "N");
    }

    // Test case 2: Another valid atom (chain A, residue 1, atom CA)
    {
        std::string cid = "//A/1/CA";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(trim_atom_name(atom->name), "CA");
    }

    // Test case 3: Invalid CID (non-existent residue)
    {
        std::string cid = "//A/9999/N";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_EQ(atom, nullptr);
    }

    // Test case 4: Chain B (water molecule)
    {
        std::string cid = "//B/1/O";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(trim_atom_name(atom->name), "O");
    }

    // Test case 5: Partial CID (chain/residue only, no atom name — returns first match)
    {
        std::string cid = "//A/10";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(trim_atom_name(atom->name), "N");
    }

    // Test case 6: Empty string CID (returns first match)
    {
        std::string cid = "";
        gemmi::Atom *atom = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(trim_atom_name(atom->name), "N");
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}