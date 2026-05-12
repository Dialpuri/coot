#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetGroup) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    {
        gemmi::Residue* found_res = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        if (res.seqid.num.value == 10) {
                            found_res = &res;
                            break;
                        }
                    }
                }
            }
        }
        
        ASSERT_NE(found_res, nullptr) << "Could not find residue //A/10";
        
        coot::protein_geometry pg;
        std::string group = pg.get_group_gemmi(*found_res);
        EXPECT_EQ(group, "peptide");
    }

    {
        gemmi::Residue* found_res = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        if (res.seqid.num.value == 999) {
                            found_res = &res;
                            break;
                        }
                    }
                }
            }
        }
        
        if (found_res) {
            coot::protein_geometry pg;
            std::string group = pg.get_group_gemmi(*found_res);
            EXPECT_FALSE(group.empty());
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}