#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <map>

TEST(OracleTest, side_chain_densities_get_block_gemmi) {
    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Collect all residues as raw pointers (stable as long as structure lives)
    std::vector<gemmi::Residue*> residues;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            residues.push_back(&res);
        }
    }

    EXPECT_FALSE(residues.empty());

    // Build a manual cache mirroring what fill_residue_blocks would produce
    std::map<gemmi::Residue*, coot::density_box_t> cache;

    // Case 1: first residue (chain A, MET, seq 1)
    {
        gemmi::Residue* res = residues[0];
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "MET");

        // Create density box using default constructor
        coot::density_box_t box;
        box.n_steps = 45;
        box.mean = 0.0;
        box.var = 0.0;
        box.is_weird = true; // all zeros => weird
        cache[res] = box;
    }

    // Case 2: middle residue (index 134)
    {
        int idx = residues.size() / 2;
        EXPECT_EQ(idx, 134);
        gemmi::Residue* res = residues[idx];
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "LEU");
        EXPECT_EQ(res->seqid.num.value, 166);

        coot::density_box_t box;
        box.n_steps = 45;
        box.mean = 0.0;
        box.var = 0.0;
        box.is_weird = true;
        cache[res] = box;
    }

    // Case 3: last residue (index 267)
    {
        int idx = residues.size() - 1;
        EXPECT_EQ(idx, 267);
        gemmi::Residue* res = residues[idx];
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "HOH");
        EXPECT_EQ(res->seqid.num.value, 1);

        coot::density_box_t box;
        box.n_steps = 45;
        box.mean = 0.0;
        box.var = 0.0;
        box.is_weird = false; // HOH with zeros => not weird (different logic in normalize)
        cache[res] = box;
    }

    // Add remaining residues with default boxes (to reach ~268 entries)
    for (gemmi::Residue* res : residues) {
        if (cache.find(res) == cache.end()) {
            coot::density_box_t box;
            box.n_steps = 45;
            box.mean = 0.0;
            box.var = 0.0;
            box.is_weird = true;
            cache[res] = box;
        }
    }

    // Check cache_size matches original
    EXPECT_EQ(cache.size(), residues.size());

    // Now test get_block_gemmi — the function under test
    // Case 1: first residue
    {
        gemmi::Residue* res = residues[0];
        coot::density_box_t block = coot::side_chain_densities::get_block_gemmi(cache, res);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "MET");
        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_NEAR(block.mean, 0.0, 1e-4);
        EXPECT_NEAR(block.var, 0.0, 1e-4);
        EXPECT_TRUE(block.is_weird);
    }

    // Case 2: middle residue
    {
        int idx = residues.size() / 2;
        EXPECT_EQ(idx, 134);
        gemmi::Residue* res = residues[idx];
        coot::density_box_t block = coot::side_chain_densities::get_block_gemmi(cache, res);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "LEU");
        EXPECT_EQ(res->seqid.num.value, 166);
        EXPECT_NEAR(block.mean, 0.0, 1e-4);
        EXPECT_NEAR(block.var, 0.0, 1e-4);
        EXPECT_TRUE(block.is_weird);
    }

    // Case 3: last residue
    {
        int idx = residues.size() - 1;
        EXPECT_EQ(idx, 267);
        gemmi::Residue* res = residues[idx];
        coot::density_box_t block = coot::side_chain_densities::get_block_gemmi(cache, res);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "HOH");
        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_NEAR(block.mean, 0.0, 1e-4);
        EXPECT_NEAR(block.var, 0.0, 1e-4);
        EXPECT_FALSE(block.is_weird);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}