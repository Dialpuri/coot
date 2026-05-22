#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(NextResidueTest, NextResidueFunction) {
    // case 1: middle residue (typical use)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1u) << "PDB loading failed";

        gemmi::CRA target_cra{nullptr, nullptr, nullptr};
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (chain.name == "A" && res.seqid.num.value == 50)
                        target_cra = gemmi::CRA{&chain, &res, nullptr};

        ASSERT_NE(target_cra.residue, nullptr) << "Could not find residue 50 in chain A";
        EXPECT_EQ(target_cra.chain->name, "A");
        EXPECT_EQ(target_cra.residue->seqid.num.value, 50);
        EXPECT_EQ(target_cra.residue->name, "ARG");

        gemmi::CRA next_cra = coot::util::next_residue_gemmi(target_cra);
        EXPECT_NE(next_cra.residue, nullptr) << "next_residue should not be nullptr for middle residue";

        EXPECT_EQ(next_cra.chain->name, "A");
        EXPECT_EQ(next_cra.residue->seqid.num.value, 51);
        EXPECT_EQ(next_cra.residue->name, "GLU");
    }

    // case 2: last residue in chain (edge case — should return nullptr)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1u) << "PDB loading failed";

        gemmi::CRA target_cra{nullptr, nullptr, nullptr};
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (chain.name == "A" && res.seqid.num.value == 298)
                        target_cra = gemmi::CRA{&chain, &res, nullptr};

        ASSERT_NE(target_cra.residue, nullptr) << "Could not find residue 298 in chain A";
        EXPECT_EQ(target_cra.chain->name, "A");
        EXPECT_EQ(target_cra.residue->seqid.num.value, 298);
        EXPECT_EQ(target_cra.residue->name, "LEU");

        gemmi::CRA next_cra = coot::util::next_residue_gemmi(target_cra);
        EXPECT_EQ(next_cra.residue, nullptr) << "next_residue should be nullptr for last residue in chain";
    }

    // case 3: null input (guard clause)
    {
        gemmi::CRA null_cra{nullptr, nullptr, nullptr};
        gemmi::CRA next_cra = coot::util::next_residue_gemmi(null_cra);
        EXPECT_EQ(next_cra.residue, nullptr) << "next_residue should be nullptr for nullptr input";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}