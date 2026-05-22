#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residues_in_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: valid chain ID "A"
    {
        std::string chain_id = "A";
        std::vector<gemmi::Residue*> v = coot::util::residues_in_chain_gemmi(&st, chain_id);
        EXPECT_EQ(v.size(), 267);
        if (!v.empty()) {
            gemmi::Residue* r = v[0];
            EXPECT_EQ(r->name, "MET");
            EXPECT_EQ(r->seqid.num.value, 1);
        }
    }

    // Case 2: valid chain ID "B"
    {
        std::string chain_id = "B";
        std::vector<gemmi::Residue*> v = coot::util::residues_in_chain_gemmi(&st, chain_id);
        EXPECT_EQ(v.size(), 1);
    }

    // Case 3: invalid chain ID "Z" (nonexistent)
    {
        std::string chain_id = "Z";
        std::vector<gemmi::Residue*> v = coot::util::residues_in_chain_gemmi(&st, chain_id);
        EXPECT_EQ(v.size(), 0);
    }

    // Case 4: null structure
    {
        std::string chain_id = "A";
        std::vector<gemmi::Residue*> v = coot::util::residues_in_chain_gemmi(nullptr, chain_id);
        EXPECT_EQ(v.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}