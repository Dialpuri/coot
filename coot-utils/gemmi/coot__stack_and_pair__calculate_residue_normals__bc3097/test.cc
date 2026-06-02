#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, calculate_residue_normals) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Find chain A
    const gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Case 1: protein residues — expected: empty map (no nucleic acid base atoms)
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        // Get residue at index 10 (0-based) from chain
        if (chain->residues.size() > 10) {
            residues_vec.push_back(std::make_pair(true, &chain->residues[10]));
        }
        if (chain->residues.size() > 20) {
            residues_vec.push_back(std::make_pair(true, &chain->residues[20]));
        }

        EXPECT_EQ(residues_vec.size(), 2u);

        std::map<const gemmi::Residue*, clipper::Coord_orth> normal_map =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normal_map.size(), 0u);
    }

    // Case 2: empty input — expected: empty map
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;

        EXPECT_EQ(residues_vec.size(), 0u);

        std::map<const gemmi::Residue*, clipper::Coord_orth> normal_map =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normal_map.size(), 0u);
    }

    // Case 3: vector containing nullptr residue — expected: empty map (nullptrs are skipped)
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        residues_vec.push_back(std::make_pair(true, nullptr));
        if (chain->residues.size() > 10) {
            residues_vec.push_back(std::make_pair(true, &chain->residues[10]));
        }

        EXPECT_EQ(residues_vec.size(), 2u);
        EXPECT_TRUE(residues_vec[0].second == nullptr);

        std::map<const gemmi::Residue*, clipper::Coord_orth> normal_map =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_EQ(normal_map.size(), 0u);
    }

    // Case 4 (complementary): Try all residues in the chain — tests non-trivial input size
    // PDB is protein-only so normal_map is still empty, but verifies no crash with large input
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        for (const gemmi::Residue& res : chain->residues) {
            residues_vec.push_back(std::make_pair(true, &res));
        }

        std::map<const gemmi::Residue*, clipper::Coord_orth> normal_map =
            coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec);

        EXPECT_NO_THROW(coot::stack_and_pair::calculate_residue_normals_gemmi(residues_vec));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}