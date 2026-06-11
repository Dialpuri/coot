#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, pucker_analysis_info_t) {
    // Setup: load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: protein residue (MET 1 in chain A) — should throw: no nucleotide base atoms
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& residue : chain.residues) {
                    if (residue.seqid.num.value == 1 && residue.name == "MET") {
                        res = &residue;
                        break;
                    }
                }
            }
        }
        std::string altconf = "";
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->name, "MET");
            EXPECT_EQ(res->seqid.num.value, 1);
        }
        EXPECT_THROW({
            coot::pucker_analysis_info_t pai;
            pai.pucker_analysis_info_t_gemmi(*res, altconf);
        }, std::runtime_error);
    }

    // Case 2: second protein residue (GLU 2 in chain A) — same exception path
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& residue : chain.residues) {
                    if (residue.seqid.num.value == 2 && residue.name == "GLU") {
                        res = &residue;
                        break;
                    }
                }
            }
        }
        std::string altconf = "";
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->name, "GLU");
            EXPECT_EQ(res->seqid.num.value, 2);
        }
        EXPECT_THROW({
            coot::pucker_analysis_info_t pai;
            pai.pucker_analysis_info_t_gemmi(*res, altconf);
        }, std::runtime_error);
    }

    // Case 3: non-existent residue — returns nullptr, skip constructor
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& residue : chain.residues) {
                    if (residue.seqid.num.value == 9999) {
                        res = &residue;
                        break;
                    }
                }
            }
        }
        std::string altconf = "";
        EXPECT_EQ(res, nullptr);
        if (res) {
            // This block should not be reached — but if it is, it should throw
            EXPECT_THROW({
                coot::pucker_analysis_info_t pai;
                pai.pucker_analysis_info_t_gemmi(*res, altconf);
            }, std::runtime_error);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
