#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <optional>
#include "function.hh"

TEST(OracleTest, get_following_residue_gemmi) {
    // Load PDB — same setup as oracle
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    // case: valid residue (chain A, seq 10) → next is GLY at seq 11
    {
        coot::residue_spec_t rs("A", 10);
        std::optional<gemmi::CRA> following =
            coot::util::get_following_residue_gemmi(rs, st);

        EXPECT_TRUE(following.has_value());
        if (following.has_value()) {
            EXPECT_NE(following->residue, nullptr);
            std::string resname = following->residue->name;
            int seqnum = following->residue->seqid.num.value;
            EXPECT_EQ(resname, "GLY");
            EXPECT_EQ(seqnum, 11);
        }
    }

    // case: last residue in chain A (seq 298) — no following → empty optional
    {
        coot::residue_spec_t rs("A", 298);
        std::optional<gemmi::CRA> following =
            coot::util::get_following_residue_gemmi(rs, st);

        EXPECT_FALSE(following.has_value());
    }

    // case: invalid chain Z — no following → empty optional
    {
        coot::residue_spec_t rs("Z", 1);
        std::optional<gemmi::CRA> following =
            coot::util::get_following_residue_gemmi(rs, st);

        EXPECT_FALSE(following.has_value());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}