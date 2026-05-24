#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(MakeCablamMarkupsTest, BasicTest) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    std::vector<std::pair<coot::residue_spec_t, double> > residues;
    coot::residue_spec_t spec1("A", 10, "");
    residues.push_back(std::make_pair(spec1, 0.5));
    
    std::vector<coot::cablam_markup_t> markups1 = coot::make_cablam_markups_gemmi(residues, st);
    EXPECT_EQ(markups1.size(), 1);
    
    if (markups1.size() > 0) {
        EXPECT_TRUE(markups1[0].cra.residue != nullptr);
        if (markups1[0].cra.residue) {
            std::string chain_id = markups1[0].cra.chain->name;
            int seq_num = markups1[0].cra.residue->seqid.num.value;
            std::string res_name = markups1[0].cra.residue->name;
            EXPECT_EQ(chain_id, "A");
            EXPECT_EQ(seq_num, 10);
            EXPECT_EQ(res_name, "ILE");
        }
    }
    
    std::vector<std::pair<coot::residue_spec_t, double> > residues2;
    coot::residue_spec_t spec2("A", 999, "");
    residues2.push_back(std::make_pair(spec2, 0.3));
    
    std::vector<coot::cablam_markup_t> markups2 = coot::make_cablam_markups_gemmi(residues2, st);
    EXPECT_EQ(markups2.size(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}