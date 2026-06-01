#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_biggest_hetgroup) {
    // Case 1: PDB with HOH hetero residue (chain B)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::CRA cra = coot::util::get_biggest_hetgroup_gemmi(&st.models[0]);
        EXPECT_NE(cra.residue, nullptr);
        ASSERT_STREQ(cra.residue->name.c_str(), "HOH");
        ASSERT_STREQ(cra.chain->name.c_str(), "B");
        EXPECT_EQ(static_cast<int>(cra.residue->seqid.num.value), 1);
        EXPECT_EQ(static_cast<int>(cra.residue->atoms.size()), 1);
    }

    // Case 2: nullptr input — verifies the guarded path
    {
        gemmi::CRA cra = coot::util::get_biggest_hetgroup_gemmi(nullptr);
        EXPECT_EQ(cra.residue, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}