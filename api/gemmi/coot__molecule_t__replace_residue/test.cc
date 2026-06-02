#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, replace_residue) {
    // Load structure using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Create a protein_geometry object
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: valid replacement — ILE to VAL
    {
        std::string cid = "//A/10";
        std::string new_type = "VAL";
        int imol_enc = 0;

        gemmi::Residue* res = coot::cid_to_residue_gemmi(cid, model);
        std::string res_before_name = res ? res->name : "nullptr";

        EXPECT_EQ(res_before_name, "ILE");

        int status = coot::replace_residue_gemmi(cid, new_type, imol_enc, geom, model);

        gemmi::Residue* res_after = coot::cid_to_residue_gemmi(cid, model);
        std::string res_after_name = res_after ? res_after->name : "nullptr";

        EXPECT_EQ(status, 1);
        EXPECT_EQ(res_after_name, "VAL");
    }

    // Case 2: invalid new residue type — should return 0
    {
        std::string cid = "//A/20";
        std::string new_type = "ZZZ";
        int imol_enc = 0;

        int status = coot::replace_residue_gemmi(cid, new_type, imol_enc, geom, model);

        EXPECT_EQ(status, 0);
    }

    // Case 3: invalid CID — residue doesn't exist
    {
        std::string cid = "//A/99999";
        std::string new_type = "GLY";
        int imol_enc = 0;

        int status = coot::replace_residue_gemmi(cid, new_type, imol_enc, geom, model);

        EXPECT_EQ(status, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}