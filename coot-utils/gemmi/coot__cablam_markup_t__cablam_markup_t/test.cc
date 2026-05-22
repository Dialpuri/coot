#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cablam_markup_t_Constructors) {
    // case 1: default constructor — no arguments
    {
        coot::cablam_markup_t_gemmi cm;
        EXPECT_EQ(cm.score, -1);
        EXPECT_EQ(cm.residue, nullptr);
    }

    // case 2: constructor with gemmi residue pointer
    {
        // Create a minimal gemmi structure for testing
        gemmi::Structure st;
        gemmi::Model model;
        gemmi::Chain chain("A");
        gemmi::Residue res(gemmi::ResidueId{gemmi::SeqId(42, ' '), "", "ALA"});
        res.name = "ALA";
        chain.residues.push_back(res);
        model.chains.push_back(chain);
        st.models.push_back(model);
        
        coot::cablam_markup_t_gemmi cm(&st.models[0].chains[0].residues[0]);
        EXPECT_EQ(cm.score, -1);
        EXPECT_EQ(cm.residue, &st.models[0].chains[0].residues[0]);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}