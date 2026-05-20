#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, PuckerAnalysisTest) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Model& model = st.models[0];
    
    // Test 1: Try to get phosphate distance on a protein residue (should fail)
    {
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr);
        
        std::string alt_conf = "";
        try {
            coot::pucker_analysis_info_t pai(res, alt_conf);
            
            gemmi::Residue* next_res = nullptr;
            for (auto& chain : model.chains) {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 11) {
                        next_res = &r;
                        break;
                    }
                }
                if (next_res) break;
            }
            
            if (next_res) {
                float dist = pai.phosphate_distance_to_base_plane_gemmi(next_res);
                FAIL() << "Expected exception for protein residue, but got distance: " << dist;
            } else {
                FAIL() << "Expected exception for protein residue, but next_res was not found";
            }
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("base lsq plane has no value"), std::string::npos);
        }
    }

    // Test 2: Try with a residue that has no following residue
    {
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.seqid.num.value == 298) {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr);
        
        std::string alt_conf = "";
        try {
            coot::pucker_analysis_info_t pai(res, alt_conf);
            gemmi::Residue* next_res = nullptr;
            float dist = pai.phosphate_distance_to_base_plane_gemmi(next_res);
            FAIL() << "Expected exception for null next_res, but got distance: " << dist;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("base lsq plane has no value"), std::string::npos);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}