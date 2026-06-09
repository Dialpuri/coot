#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <map>
#include <vector>
#include <string>
#include <utility>

static std::string trim_atom_name(const std::string& s) {
    std::string trimmed = s;
    size_t start = trimmed.find_first_not_of(" ");
    if (start != std::string::npos) {
        trimmed = trimmed.substr(start);
    } else {
        return "";
    }
    size_t end = trimmed.find_last_not_of(" ");
    if (end != std::string::npos) {
        trimmed = trimmed.substr(0, end + 1);
    }
    return trimmed;
}

TEST(OracleTest, make_neighbour_map) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Helper: collect CRAs for atoms in a given chain/residue
    auto collect_cras = [&](const std::string& chain_id, int seqnum) {
        std::vector<gemmi::CRA> cras;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != seqnum) continue;
                for (gemmi::Atom& atom : res.atoms) {
                    cras.push_back(gemmi::CRA{&chain, &res, &atom});
                }
            }
        }
        return cras;
    };

    // Case 1: Selection covering all atoms in chain A residue 10 (ILE)
    {
        auto cras_1 = collect_cras("A", 10);
        auto cras_2 = cras_1;

        EXPECT_EQ(cras_1.size(), 8u);

        std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>> neighbour_map =
            coot::h_bonds::make_neighbour_map_gemmi(cras_1, cras_2);

        EXPECT_EQ(neighbour_map.size(), 8u);

        bool found_N = false;
        bool found_CA = false;
        bool found_C = false;
        bool found_O = false;
        bool found_CB = false;
        bool found_CG1 = false;
        bool found_CG2 = false;
        bool found_CD1 = false;

        for (auto const &entry : neighbour_map) {
            gemmi::Atom *at = entry.first;
            std::string atom_name = trim_atom_name(at->name);

            // Find the CRA for this atom to verify residue info
            gemmi::Residue* res = nullptr;
            gemmi::Chain* ch = nullptr;
            for (auto const& cra : cras_1) {
                if (cra.atom == at) {
                    res = cra.residue;
                    ch = cra.chain;
                    break;
                }
            }
            ASSERT_NE(res, nullptr);
            ASSERT_NE(ch, nullptr);

            // Verify residue info
            EXPECT_STREQ(ch->name.c_str(), "A");
            EXPECT_EQ(res->seqid.num.value, 10);
            EXPECT_STREQ(res->name.c_str(), "ILE");

            if (atom_name == "N") {
                found_N = true;
                EXPECT_EQ(entry.second.size(), 1u);
                EXPECT_NEAR(entry.second[0].second, 1.45564f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "CA");
            } else if (atom_name == "CA") {
                found_CA = true;
                EXPECT_EQ(entry.second.size(), 3u);
                EXPECT_NEAR(entry.second[0].second, 1.45564f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "N");
                EXPECT_NEAR(entry.second[1].second, 1.52938f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[1].first->name), "C");
                EXPECT_NEAR(entry.second[2].second, 1.54292f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[2].first->name), "CB");
            } else if (atom_name == "C") {
                found_C = true;
                EXPECT_EQ(entry.second.size(), 2u);
                EXPECT_NEAR(entry.second[0].second, 1.23193f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "O");
                EXPECT_NEAR(entry.second[1].second, 1.52938f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[1].first->name), "CA");
            } else if (atom_name == "O") {
                found_O = true;
                EXPECT_EQ(entry.second.size(), 1u);
                EXPECT_NEAR(entry.second[0].second, 1.23193f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "C");
            } else if (atom_name == "CB") {
                found_CB = true;
                EXPECT_EQ(entry.second.size(), 3u);
                EXPECT_NEAR(entry.second[0].second, 1.52723f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "CG2");
                EXPECT_NEAR(entry.second[1].second, 1.53186f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[1].first->name), "CG1");
                EXPECT_NEAR(entry.second[2].second, 1.54292f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[2].first->name), "CA");
            } else if (atom_name == "CG1") {
                found_CG1 = true;
                EXPECT_EQ(entry.second.size(), 2u);
                EXPECT_NEAR(entry.second[0].second, 1.52069f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "CD1");
                EXPECT_NEAR(entry.second[1].second, 1.53186f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[1].first->name), "CB");
            } else if (atom_name == "CG2") {
                found_CG2 = true;
                EXPECT_EQ(entry.second.size(), 1u);
                EXPECT_NEAR(entry.second[0].second, 1.52723f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "CB");
            } else if (atom_name == "CD1") {
                found_CD1 = true;
                EXPECT_EQ(entry.second.size(), 1u);
                EXPECT_NEAR(entry.second[0].second, 1.52069f, 1e-4);
                EXPECT_EQ(trim_atom_name(entry.second[0].first->name), "CG1");
            }
        }

        EXPECT_TRUE(found_N);
        EXPECT_TRUE(found_CA);
        EXPECT_TRUE(found_C);
        EXPECT_TRUE(found_O);
        EXPECT_TRUE(found_CB);
        EXPECT_TRUE(found_CG1);
        EXPECT_TRUE(found_CG2);
        EXPECT_TRUE(found_CD1);
    }

    // Case 2: Empty selection (non-existent residue) — should produce empty map
    {
        std::vector<gemmi::CRA> cras_1; // empty
        std::vector<gemmi::CRA> cras_2; // empty

        EXPECT_EQ(cras_1.size(), 0u);

        std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>> neighbour_map =
            coot::h_bonds::make_neighbour_map_gemmi(cras_1, cras_2);

        EXPECT_EQ(neighbour_map.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
