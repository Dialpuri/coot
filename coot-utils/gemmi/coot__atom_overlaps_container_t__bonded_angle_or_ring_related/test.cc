#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <iostream>
#include <string>
#include <algorithm>

static std::string interaction_name(coot::bonded_atom_interaction_type ait) {
    switch (ait) {
        case coot::BAI_CLASHABLE: return "CLASHABLE";
        case coot::BAI_BONDED:    return "BONDED";
        case coot::BAI_IGNORED:   return "IGNORED";
    }
    return "UNKNOWN";
}

TEST(OracleTest, bonded_angle_or_ring_related) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model &model = st.models[0];

    // Find atoms by name: N/CA/C/N2/CA2/CB2/CB1/CB3
    gemmi::CRA cra_N1{nullptr, nullptr, nullptr};  // A/1/N
    gemmi::CRA cra_CA1{nullptr, nullptr, nullptr}; // A/1/CA
    gemmi::CRA cra_C1{nullptr, nullptr, nullptr};  // A/1/C
    gemmi::CRA cra_N2{nullptr, nullptr, nullptr};  // A/2/N
    gemmi::CRA cra_CA2{nullptr, nullptr, nullptr}; // A/2/CA
    gemmi::CRA cra_CB2{nullptr, nullptr, nullptr}; // A/2/CB
    gemmi::CRA cra_CB1{nullptr, nullptr, nullptr}; // A/1/CB (sidechain)
    gemmi::CRA cra_CB3{nullptr, nullptr, nullptr}; // A/3/CB (sidechain, non-adjacent test)

    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue &res : chain.residues) {
            int seq = res.seqid.num.value;
            for (gemmi::Atom &atom : res.atoms) {
                std::string aname = atom.name;

                // Trim trailing AND leading spaces for comparison
                std::string trimmed = aname;
                size_t first = trimmed.find_first_not_of(" ");
                if (first != std::string::npos) trimmed = trimmed.substr(first);
                size_t last = trimmed.find_last_not_of(" ");
                if (last != std::string::npos) trimmed = trimmed.substr(0, last + 1);

                gemmi::CRA cra{&chain, &res, &atom};
                if (seq == 1) {
                    if (trimmed == "N")  cra_N1 = cra;
                    if (trimmed == "CA") cra_CA1 = cra;
                    if (trimmed == "C")  cra_C1 = cra;
                    if (trimmed == "CB") cra_CB1 = cra;
                }
                if (seq == 2) {
                    if (trimmed == "N")  cra_N2 = cra;
                    if (trimmed == "CA") cra_CA2 = cra;
                    if (trimmed == "CB") cra_CB2 = cra;
                }
                if (seq == 3) {
                    if (trimmed == "CB") cra_CB3 = cra;
                }
            }
        }
    }

    // Get connections from structure
    const std::vector<gemmi::Connection> &connections = st.connections;

    // Case 1: Same residue bonded atoms (N-CA in ALA 1)
    {
        ASSERT_NE(cra_N1.atom, nullptr);
        ASSERT_NE(cra_CA1.atom, nullptr);
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> bonded_neighbours;
        std::map<std::string, std::vector<std::vector<std::string>>> ring_list_map;

        auto ait = coot::bonded_angle_or_ring_related_gemmi(cra_N1, cra_CA1, false,
                                                            &bonded_neighbours, &ring_list_map,
                                                            &model, connections,
                                                            &geom, false);

        EXPECT_EQ(interaction_name(ait), "BONDED");
    }

    // Case 2: Adjacent residue bonded atoms (C of res1, N of res2)
    {
        ASSERT_NE(cra_C1.atom, nullptr);
        ASSERT_NE(cra_N2.atom, nullptr);
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> bonded_neighbours;
        std::map<std::string, std::vector<std::vector<std::string>>> ring_list_map;

        auto ait = coot::bonded_angle_or_ring_related_gemmi(cra_C1, cra_N2, false,
                                                            &bonded_neighbours, &ring_list_map,
                                                            &model, connections,
                                                            &geom, false);

        EXPECT_EQ(interaction_name(ait), "BONDED");
    }

    // Case 3: Same residue 1-3 related atoms (CA-CB in ALA 2)
    {
        ASSERT_NE(cra_CA2.atom, nullptr);
        ASSERT_NE(cra_CB2.atom, nullptr);
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> bonded_neighbours;
        std::map<std::string, std::vector<std::vector<std::string>>> ring_list_map;

        auto ait = coot::bonded_angle_or_ring_related_gemmi(cra_CA2, cra_CB2, false,
                                                            &bonded_neighbours, &ring_list_map,
                                                            &model, connections,
                                                            &geom, false);

        EXPECT_EQ(interaction_name(ait), "BONDED");
    }

    // Case 4: Non-adjacent residues (should be CLASHABLE)
    {
        ASSERT_NE(cra_CB1.atom, nullptr);
        ASSERT_NE(cra_CB3.atom, nullptr);
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> bonded_neighbours;
        std::map<std::string, std::vector<std::vector<std::string>>> ring_list_map;

        auto ait = coot::bonded_angle_or_ring_related_gemmi(cra_CB1, cra_CB3, false,
                                                            &bonded_neighbours, &ring_list_map,
                                                            &model, connections,
                                                            &geom, false);

        EXPECT_EQ(interaction_name(ait), "CLASHABLE");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
