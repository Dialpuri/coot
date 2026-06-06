#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number in the first model
static gemmi::Residue* find_residue(gemmi::Model& model, std::string chain_id, int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, add_tetrahedral_hydrogen) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: THR HB on CB with neighbours CA, OG1, CG2
    {
        gemmi::Residue* res = find_residue(model, "A", 14);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "THR");

        int atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_before, 14);

        std::vector<std::string> neighbours = {"CA", "OG1", "CG2"};
        coot::reduce::add_tetrahedral_hydrogen_gemmi("HB", "CB", neighbours, 1.0, *res);

        int atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_after, 15);
        EXPECT_EQ(atoms_after - atoms_before, 1);
    }

    // Case 2: VAL HB on CB with neighbours CA, CG1, CG2
    {
        gemmi::Residue* res = find_residue(model, "A", 7);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "VAL");

        int atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_before, 16);

        std::vector<std::string> neighbours = {"CA", "CG1", "CG2"};
        coot::reduce::add_tetrahedral_hydrogen_gemmi("HB", "CB", neighbours, 1.0, *res);

        int atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_after, 17);
        EXPECT_EQ(atoms_after - atoms_before, 1);
    }

    // Case 3: Edge - wrong number of neighbours (not 3) -> WARNING path
    {
        gemmi::Residue* res = find_residue(model, "A", 5);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "GLN");

        int atoms_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_before, 17);

        std::vector<std::string> neighbours = {"CA", "CG1"};
        coot::reduce::add_tetrahedral_hydrogen_gemmi("HB", "CB", neighbours, 1.0, *res);

        int atoms_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(atoms_after, 17);
        EXPECT_EQ(atoms_after - atoms_before, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}