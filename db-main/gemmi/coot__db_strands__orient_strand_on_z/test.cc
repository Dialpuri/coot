#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, orient_strand_on_z) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: select a few residues in chain A (residues 8-14)
    {
        int SelHnd = 1;  // Simulate selection handle
        std::vector<const gemmi::Residue*> SelResidues;

        // Find chain A and collect residues 8-14
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq >= 8 && seq <= 14) {
                            SelResidues.push_back(&res);
                        }
                    }
                }
            }
        }

        int nSelResidues = static_cast<int>(SelResidues.size());

        EXPECT_EQ(SelHnd, 1);
        EXPECT_EQ(nSelResidues, 7);

        std::pair<bool, clipper::RTop_orth> result = coot::db_strands::orient_strand_on_z_gemmi(SelResidues);

        EXPECT_TRUE(result.first);
    }

    // Case 2: select a single residue in chain A (residue 10)
    {
        int SelHnd = 1;
        std::vector<const gemmi::Residue*> SelResidues;

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq == 10) {
                            SelResidues.push_back(&res);
                        }
                    }
                }
            }
        }

        int nSelResidues = static_cast<int>(SelResidues.size());

        EXPECT_EQ(SelHnd, 1);
        EXPECT_EQ(nSelResidues, 1);

        std::pair<bool, clipper::RTop_orth> result = coot::db_strands::orient_strand_on_z_gemmi(SelResidues);

        EXPECT_TRUE(result.first);
    }

    // Case 3: empty selection (non-existent residue 9999)
    {
        int SelHnd = 1;
        std::vector<const gemmi::Residue*> SelResidues;

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq == 9999) {
                            SelResidues.push_back(&res);
                        }
                    }
                }
            }
        }

        int nSelResidues = static_cast<int>(SelResidues.size());

        EXPECT_EQ(SelHnd, 1);
        EXPECT_EQ(nSelResidues, 0);

        std::pair<bool, clipper::RTop_orth> result = coot::db_strands::orient_strand_on_z_gemmi(SelResidues);

        EXPECT_TRUE(result.first);
    }

    // Case 4: HOH water residue (complementary — should fail due to missing backbone atoms)
    {
        int SelHnd = 1;
        std::vector<const gemmi::Residue*> SelResidues;

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "B") {
                    for (gemmi::Residue& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq == 1) {
                            SelResidues.push_back(&res);
                        }
                    }
                }
            }
        }

        int nSelResidues = static_cast<int>(SelResidues.size());

        EXPECT_EQ(SelHnd, 1);
        EXPECT_EQ(nSelResidues, 1);

        std::pair<bool, clipper::RTop_orth> result = coot::db_strands::orient_strand_on_z_gemmi(SelResidues);

        EXPECT_FALSE(result.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}