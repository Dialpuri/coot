#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps_container_t) {
    // Initialize protein geometry (required for any geometry-dependent work)
    coot::protein_geometry geom;
    geom.init_standard();

    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    gemmi::Model const* model = &st.models[0];

    // Helper to find a residue by chain/resseq (matching MMDB get_residue("//A/10"))
    auto find_residue = [&](std::string chain_id, int resseq) -> gemmi::Residue const* {
        for (auto& ch : model->chains) {
            if (ch.name != chain_id) continue;
            for (auto& res : ch.residues) {
                if (res.seqid.num.value == resseq) {
                    return &res;
                }
            }
        }
        return nullptr;
    };

    // Case 1: Two different neighbouring residues (typical usage)
    {
        gemmi::Residue const* res_central = find_residue("A", 10);
        gemmi::Residue const* neighbour   = find_residue("A", 11);

        EXPECT_NE(res_central, nullptr);
        EXPECT_NE(neighbour, nullptr);

        if (res_central && neighbour) {
            auto aoc = coot::atom_overlaps_container_t_gemmi(
                res_central, neighbour, model, &geom);
            EXPECT_EQ(aoc.neighbours.size(), 1u);
            EXPECT_NEAR(aoc.score(), 0.0f, 1e-4);
        }
    }

    // Case 2: Residue with itself as neighbour (edge: same residue passed twice)
    {
        gemmi::Residue const* res_central = find_residue("A", 10);

        EXPECT_NE(res_central, nullptr);

        if (res_central) {
            auto aoc = coot::atom_overlaps_container_t_gemmi(
                res_central, res_central, model, &geom);
            EXPECT_EQ(aoc.neighbours.size(), 1u);
            EXPECT_NEAR(aoc.score(), 0.0f, 1e-4);
        }
    }

    // Case 3: Two residues far apart (ALA 10 and TYR 94 — should still work, just no close contacts)
    {
        gemmi::Residue const* res_central = find_residue("A", 10);
        gemmi::Residue const* neighbour   = find_residue("A", 94);

        EXPECT_NE(res_central, nullptr);
        EXPECT_NE(neighbour, nullptr);

        if (res_central && neighbour) {
            auto aoc = coot::atom_overlaps_container_t_gemmi(
                res_central, neighbour, model, &geom);
            EXPECT_EQ(aoc.neighbours.size(), 1u);
            EXPECT_NEAR(aoc.score(), 0.0f, 1e-4);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
