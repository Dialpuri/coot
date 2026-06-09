#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a residue by chain name and sequence number in a gemmi Structure
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, chi_angles_change_by) {
    coot::protein_geometry geom;
    geom.init_standard();
    
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];
    
    // Case 1: LYS residue (has chi1) — valid chi angle
    {
        gemmi::Residue* res = find_residue(model, "A", 6);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "LYS");
        EXPECT_EQ(res->seqid.num.value, 6);

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(1, 90.0, res, &geom);
        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 0.312303, 1e-4);
    }

    // Case 2: GLU residue (has chi1, chi2, chi3) — chi3
    {
        gemmi::Residue* res = find_residue(model, "A", 8);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->name, "GLU");
        EXPECT_EQ(res->seqid.num.value, 8);

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(3, 180.0, res, &geom);
        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 1.04307, 1e-4);
    }

    // Case 3: Invalid ichi (5) — should fail with status 2
    {
        gemmi::Residue* res = find_residue(model, "A", 6);
        EXPECT_NE(res, nullptr);

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(5, 90.0, res, &geom);
        EXPECT_EQ(result.first, 2);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
