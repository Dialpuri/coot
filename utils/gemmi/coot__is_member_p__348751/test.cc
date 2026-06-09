#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <string>

TEST(OracleTest, is_member_p) {
    // Case 1: string version — member found
    {
        std::vector<std::string> v = {"ALA", "GLY", "VAL"};
        std::string a = "GLY";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_TRUE(result);
    }

    // Case 2: string version — not found
    {
        std::vector<std::string> v = {"ALA", "GLY", "VAL"};
        std::string a = "PHE";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_FALSE(result);
    }

    // Case 3: int version — member found
    {
        std::vector<int> v = {1, 5, 10, 20};
        int a = 10;
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_TRUE(result);
    }

    // Case 4: int version — not found
    {
        std::vector<int> v = {1, 5, 10, 20};
        int a = 99;
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_FALSE(result);
    }

    // Case 5: Residue* version
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        if (st.models.empty()) {
            ADD_FAILURE() << "PDB load failed";
            return;
        }

        gemmi::Model &model = st.models[0];
        if (model.chains.empty()) {
            ADD_FAILURE() << "No chains found";
            return;
        }

        gemmi::Chain &chain = model.chains[0];
        if (chain.residues.size() < 3) {
            ADD_FAILURE() << "Not enough residues";
            return;
        }

        gemmi::Residue *r1 = &(chain.residues[0]);
        gemmi::Residue *r2 = &(chain.residues[1]);
        gemmi::Residue *r3 = &(chain.residues[2]);

        std::vector<gemmi::Residue *> resvec = {r1, r2};

        // Case 5a: residue is in vector
        {
            bool result = coot::is_member_p_gemmi(resvec, r1);
            EXPECT_TRUE(result);
        }

        // Case 5b: residue not in vector
        {
            bool result = coot::is_member_p_gemmi(resvec, r3);
            EXPECT_FALSE(result);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
