#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_het_residue) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Standard amino acid residue (ILE at //A/10) - should NOT be het
    {
        bool found = false;
        for (gemmi::Model& model : st.models)
            for (gemmi::Chain& chain : model.chains)
                for (gemmi::Residue& res : chain.residues)
                    if (res.seqid.num.value == 10 && chain.name == "A") {
                        found = true;
                        bool result = coot::cho::is_het_residue_gemmi(&res);
                        EXPECT_FALSE(result) << "ILE at A/10 should not be het residue";
                        break;
                    }
        EXPECT_TRUE(found) << "Residue at A/10 not found";
    }

    // Case 2: Water residue (HOH at //B/1) - should be het
    {
        bool found = false;
        for (gemmi::Model& model : st.models)
            for (gemmi::Chain& chain : model.chains)
                for (gemmi::Residue& res : chain.residues)
                    if (res.seqid.num.value == 1 && chain.name == "B") {
                        found = true;
                        bool result = coot::cho::is_het_residue_gemmi(&res);
                        EXPECT_TRUE(result) << "HOH at B/1 should be het residue";
                        break;
                    }
        EXPECT_TRUE(found) << "Residue at B/1 not found";
    }

    // Case 3: NULL residue (edge case) - should return false
    {
        const gemmi::Residue* null_res = nullptr;
        bool result = coot::cho::is_het_residue_gemmi(null_res);
        EXPECT_FALSE(result) << "nullptr should return false";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}