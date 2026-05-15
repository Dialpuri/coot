#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(StackAndPairTest, get_base_atom_names) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: ILE residue A/10
    {
        gemmi::Residue *res = nullptr;
        std::string chain_id;
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                for (auto &r : chain.residues) {
                    if (r.name == "ILE" && r.seqid.num.value == 10 && r.seqid.icode == ' ') {
                        res = &r;
                        chain_id = chain.name;
                        break;
                    }
                }
                if (res) break;
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(chain_id.c_str(), "A");
        EXPECT_EQ(res->seqid.num.value, 10);
        std::vector<std::string> names = coot::stack_and_pair::get_base_atom_names_gemmi(*res);
        EXPECT_EQ(names.size(), 0u);
    }

    // case 2: GLY residue A/11
    {
        gemmi::Residue *res = nullptr;
        std::string chain_id;
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                for (auto &r : chain.residues) {
                    if (r.name == "GLY" && r.seqid.num.value == 11 && r.seqid.icode == ' ') {
                        res = &r;
                        chain_id = chain.name;
                        break;
                    }
                }
                if (res) break;
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(chain_id.c_str(), "A");
        EXPECT_EQ(res->seqid.num.value, 11);
        std::vector<std::string> names = coot::stack_and_pair::get_base_atom_names_gemmi(*res);
        EXPECT_EQ(names.size(), 0u);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}