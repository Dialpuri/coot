#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, PuckerAnalysisNonNucleotide) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test 1: Non-nucleotide residue ILE at position 10
    {
        gemmi::Model& model = st.models[0];
        bool found_ile = false;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    EXPECT_EQ(res.name, "ILE");
                    try {
                        coot::pucker_analysis_info_t pai(&res, "");
                        FAIL() << "Expected exception not thrown";
                    } catch (const std::exception& e) {
                        EXPECT_STREQ(e.what(), "base lsq plane has no value");
                    }
                    found_ile = true;
                    break;
                }
            }
            if (found_ile) break;
        }
    }

    // Test 2: Non-nucleotide residue LYS at position 20
    {
        gemmi::Model& model = st.models[0];
        bool found_lys = false;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 20) {
                    EXPECT_EQ(res.name, "LYS");
                    try {
                        coot::pucker_analysis_info_t pai(&res, "");
                        FAIL() << "Expected exception not thrown";
                    } catch (const std::exception& e) {
                        EXPECT_STREQ(e.what(), "base lsq plane has no value");
                    }
                    found_lys = true;
                    break;
                }
            }
            if (found_lys) break;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}