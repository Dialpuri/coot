#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_member_p) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];

    gemmi::Chain *chainA = nullptr;
    for (gemmi::Chain &ch : model.chains) {
        if (ch.name == "A") {
            chainA = &ch;
            break;
        }
    }
    ASSERT_NE(chainA, nullptr);

    // Collect first 3 residues into a vector
    std::vector<gemmi::Residue*> res_vec;
    int nres = static_cast<int>(chainA->residues.size());
    int take = (nres > 3) ? 3 : nres;
    for (int i = 0; i < take; i++) {
        gemmi::Residue *r = &(chainA->residues[static_cast<size_t>(i)]);
        res_vec.push_back(r);
    }

    // Case 1: residue that IS in the vector (ground truth: true)
    {
        gemmi::Residue *target = res_vec[0];
        bool result = coot::is_member_p_gemmi(res_vec, target);
        EXPECT_EQ(res_vec.size(), 3u);
        EXPECT_TRUE(result);
    }

    // Case 2: residue that is NOT in the vector (ground truth: false)
    {
        gemmi::Residue *target = &(chainA->residues[static_cast<size_t>(nres - 1)]);
        bool result = coot::is_member_p_gemmi(res_vec, target);
        EXPECT_EQ(res_vec.size(), 3u);
        EXPECT_FALSE(result);
    }

    // Case 3: empty vector (ground truth: false)
    {
        std::vector<gemmi::Residue*> empty_vec;
        gemmi::Residue *target = res_vec[0];
        bool result = coot::is_member_p_gemmi(empty_vec, target);
        EXPECT_EQ(empty_vec.size(), 0u);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
