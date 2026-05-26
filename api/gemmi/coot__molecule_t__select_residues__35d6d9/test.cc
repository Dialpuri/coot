#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SelectResidues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: valid range in chain A (should return multiple residues)
    {
        std::string chain_id = "A";
        int resno_start = 10;
        int resno_end = 15;
        std::vector<gemmi::CRA> rv = coot::molecule_t::select_residues_gemmi(chain_id, resno_start, resno_end, st);
        EXPECT_EQ(rv.size(), 6u);
        if (!rv.empty()) {
            EXPECT_EQ(rv[0].residue->name, "ILE");
            EXPECT_EQ(rv[0].chain->name, "A");
            EXPECT_EQ(rv[0].residue->seqid.num.value, 10);
            EXPECT_EQ(rv.back().residue->name, "TYR");
            EXPECT_EQ(rv.back().chain->name, "A");
            EXPECT_EQ(rv.back().residue->seqid.num.value, 15);
        }
    }

    // case 2: single residue range (edge case)
    {
        std::string chain_id = "A";
        int resno_start = 50;
        int resno_end = 50;
        std::vector<gemmi::CRA> rv = coot::molecule_t::select_residues_gemmi(chain_id, resno_start, resno_end, st);
        EXPECT_EQ(rv.size(), 1u);
        if (!rv.empty()) {
            EXPECT_EQ(rv[0].residue->name, "ARG");
            EXPECT_EQ(rv[0].chain->name, "A");
            EXPECT_EQ(rv[0].residue->seqid.num.value, 50);
        }
    }

    // case 3: empty range (resno_start > resno_end)
    {
        std::string chain_id = "A";
        int resno_start = 100;
        int resno_end = 90;
        std::vector<gemmi::CRA> rv = coot::molecule_t::select_residues_gemmi(chain_id, resno_start, resno_end, st);
        EXPECT_EQ(rv.size(), 0u);
    }

    // case 4: non-existent chain
    {
        std::string chain_id = "Z";
        int resno_start = 10;
        int resno_end = 20;
        std::vector<gemmi::CRA> rv = coot::molecule_t::select_residues_gemmi(chain_id, resno_start, resno_end, st);
        EXPECT_EQ(rv.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}