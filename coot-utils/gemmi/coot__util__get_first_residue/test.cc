#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetFirstResidue) {
    // case 1: valid molecule — should find first residue
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::CRA cra = coot::util::get_first_residue_gemmi(st);
        EXPECT_NE(cra.residue, nullptr);

        const std::string& chain_id = cra.chain->name;
        int seqnum = cra.residue->seqid.num.value;
        const std::string& resname = cra.residue->name;
        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(resname, "MET");
    }

    // case 2: empty structure — no model
    {
        gemmi::Structure st;
        gemmi::CRA cra = coot::util::get_first_residue_gemmi(st);
        EXPECT_EQ(cra.residue, nullptr);
    }

    // case 3: empty model — no chains
    {
        gemmi::Structure st;
        gemmi::Model model;
        st.models.push_back(model);
        gemmi::CRA cra = coot::util::get_first_residue_gemmi(st);
        EXPECT_EQ(cra.residue, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}