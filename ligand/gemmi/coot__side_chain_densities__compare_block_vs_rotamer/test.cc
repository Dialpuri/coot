#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SideChainDensitiesTest) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.name == "TYR" && r.seqid.num.value == 77 && chain.name == "A") {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        if (res) break;
    }
    
    ASSERT_NE(res, nullptr) << "Residue not found";
    
    EXPECT_EQ(res->seqid.num.value, 77);
    EXPECT_EQ(res->name, "TYR");
    
    coot::density_box_t db;
    db.residue_p = res;
    
    coot::density_box_t empty_db;
    empty_db.residue_p = nullptr;
    
    clipper::Xmap<float> default_xmap;
    std::string rotamer_dir = "/nonexistent/path";
    
    std::pair<bool, double> result = coot::compare_block_vs_rotamer_gemmi(empty_db, res, rotamer_dir, default_xmap);
    
    EXPECT_EQ(result.first, false);
    EXPECT_NEAR(result.second, 0.0, 1e-10);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}