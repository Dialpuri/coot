#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetResidueAndCAPosition) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue with CA atom (first residue in chain A)
    {
        bool found = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 1) {
                            found = true;
                            auto result = coot::util::get_CA_position_in_residue_gemmi(res);
                            EXPECT_TRUE(result.first);
                            EXPECT_NEAR(result.second.x(), 13.401, 1e-3);
                            EXPECT_NEAR(result.second.y(), 3.78, 1e-3);
                            EXPECT_NEAR(result.second.z(), 74.703, 1e-3);
                            
                            EXPECT_EQ(chain.name, "A");
                            EXPECT_EQ(res.seqid.num.value, 1);
                            EXPECT_EQ(res.name, "MET");
                        }
        EXPECT_TRUE(found) << "Residue not found";
    }

    // case 2: another valid residue (residue 10 in chain A)
    {
        bool found = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10) {
                            found = true;
                            auto result = coot::util::get_CA_position_in_residue_gemmi(res);
                            EXPECT_TRUE(result.first);
                            EXPECT_NEAR(result.second.x(), 27.175, 1e-3);
                            EXPECT_NEAR(result.second.y(), 12.237, 1e-3);
                            EXPECT_NEAR(result.second.z(), 67.71, 1e-3);
                            
                            EXPECT_EQ(chain.name, "A");
                            EXPECT_EQ(res.seqid.num.value, 10);
                            EXPECT_EQ(res.name, "ILE");
                        }
        EXPECT_TRUE(found) << "Residue not found";
    }

    // case 3: residue with different chain (chain B) - should not find residue 5 (water only)
    {
        bool found = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "B")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 5)
                            found = true;
        EXPECT_FALSE(found) << "Residue should not be found for chain B / 5";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}