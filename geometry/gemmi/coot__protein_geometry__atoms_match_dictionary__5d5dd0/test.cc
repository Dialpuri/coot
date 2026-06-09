#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atoms_match_dictionary) {
    // Read PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Find chain A, residue ALA 21
    gemmi::Model& model = st.models[0];
    gemmi::Residue* res_p = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.name == "ALA" && res.seqid.num.value == 21) {
                    res_p = &res;
                    break;
                }
            }
        }
    }

    // Case 1: Standard residue (ALA 21), check_hydrogens_too_flag=false
    {
        ASSERT_NE(res_p, nullptr);
        EXPECT_EQ(std::string(res_p->name), "ALA");

        auto result = coot::atoms_match_dictionary_gemmi(*res_p, false, false);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 5u);
        EXPECT_EQ(result.second[0], " N  ");
        EXPECT_EQ(result.second[1], " CA ");
        EXPECT_EQ(result.second[2], " C  ");
        EXPECT_EQ(result.second[3], " O  ");
        EXPECT_EQ(result.second[4], " CB ");
    }

    // Case 2: Same residue but checking hydrogens
    {
        ASSERT_NE(res_p, nullptr);

        auto result = coot::atoms_match_dictionary_gemmi(*res_p, true, false);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 10u);
        EXPECT_EQ(result.second[0], " N  ");
        EXPECT_EQ(result.second[1], " CA ");
        EXPECT_EQ(result.second[2], " C  ");
        EXPECT_EQ(result.second[3], " O  ");
        EXPECT_EQ(result.second[4], " CB ");
        EXPECT_EQ(result.second[5], " H  ");
        EXPECT_EQ(result.second[6], " HA ");
        EXPECT_EQ(result.second[7], " HB3");
        EXPECT_EQ(result.second[8], " HB2");
        EXPECT_EQ(result.second[9], " HB1");
    }

    // Case 3: Non-existent residue (sequence number 9999)
    {
        gemmi::Residue* res9999 = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 9999) {
                        res9999 = &res;
                        break;
                    }
                }
            }
        }
        EXPECT_EQ(res9999, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}