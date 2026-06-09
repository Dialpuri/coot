#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_group) {
    // Initialize protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: standard residue ILE (oracle prints res_name: ILE)
    {
        // Find residue ILE at chain A, seqnum 10
        const gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : st.models[0].chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 10 && r.name == "ILE") {
                        res = &r;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "ILE");
        std::string group = coot::get_group_gemmi(geom, *res);
        EXPECT_EQ(group, "peptide");
    }

    // Case 2: standard residue GLY (oracle prints res_name: GLY)
    {
        // Find residue GLY at chain A, seqnum 11
        const gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : st.models[0].chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 11 && r.name == "GLY") {
                        res = &r;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "GLY");
        std::string group = coot::get_group_gemmi(geom, *res);
        EXPECT_EQ(group, "peptide");
    }

    // Case 3: invalid residue (non-existent)
    {
        const gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : st.models[0].chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 9999) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        EXPECT_EQ(res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
