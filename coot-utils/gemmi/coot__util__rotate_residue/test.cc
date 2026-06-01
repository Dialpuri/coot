#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

// Simple CID parser that validates chain ID, sequence number, AND residue name
// Returns nullptr when any component doesn't match (like the MMDB cid_to_residue does).
static gemmi::Residue* find_residue_by_cid(gemmi::Model& model, const std::string& cid) {
    // Parse CID like "//A/10/ALA" → chain "A", seqnum 10, resname "ALA"
    std::istringstream iss(cid);
    std::string token;
    std::vector<std::string> tokens;

    // Split by '/'
    while (std::getline(iss, token, '/')) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }

    if (tokens.size() != 3) return nullptr;

    std::string chain_id_str = tokens[0];
    int seqnum = std::stoi(tokens[1]);
    std::string resname = tokens[2];

    char chain_id = chain_id_str[0];

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == std::string(1, chain_id)) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seqnum && res.name == resname) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, rotate_residue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Valid residue path — but cid uses wrong residue name, so cid_to_residue returns nullptr
    {
        std::string cid = "//A/10/ALA";
        gemmi::Residue* residue_p = find_residue_by_cid(model, cid);
        EXPECT_EQ(residue_p, nullptr);
    }

    // Case 2: nullptr residue — calling rotate_residue with nullptr should not crash
    {
        std::string cid = "//Z/9999/XXX";
        gemmi::Residue* residue_p = find_residue_by_cid(model, cid);
        EXPECT_EQ(residue_p, nullptr);

        clipper::Coord_orth direction(0, 0, 1);
        clipper::Coord_orth origin_shift(0, 0, 0);
        EXPECT_NO_THROW(coot::util::rotate_residue_gemmi(residue_p, direction, origin_shift, M_PI));
    }

    // Case 3: Zero angle — but cid uses wrong residue name, so cid_to_residue returns nullptr
    {
        std::string cid = "//A/11/ARG";
        gemmi::Residue* residue_p = find_residue_by_cid(model, cid);
        EXPECT_EQ(residue_p, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}