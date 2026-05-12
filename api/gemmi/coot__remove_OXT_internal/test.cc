#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test fixture path (use VERBATIM from requirements)
static const std::string TEST_PDB_PATH = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, remove_OXT_internal) {
    // Load the structure
    gemmi::Structure st = gemmi::read_pdb_file(TEST_PDB_PATH);
    ASSERT_FALSE(st.models.empty());

    // Case 1: Terminal C-terminus residue (should have OXT, which gets removed)
    {
        // Target: chain A, residue 298 (last residue in chain A)
        gemmi::Residue* target_res = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 298)
                            target_res = &res;

        ASSERT_NE(target_res, nullptr) << "residue //A/298 not found";

        int oxt_count_before = 0;
        for (const auto& atom : target_res->atoms) {
            std::string name = atom.name;
            while (!name.empty() && name.back() == ' ')
                name.pop_back();
            if (name == "OXT")
                oxt_count_before++;
        }

        EXPECT_EQ(oxt_count_before, 1);

        // Call the function
        coot::remove_OXT_internal_gemmi(target_res);

        int oxt_count_after = 0;
        for (const auto& atom : target_res->atoms) {
            std::string name = atom.name;
            while (!name.empty() && name.back() == ' ')
                name.pop_back();
            if (name == "OXT")
                oxt_count_after++;
        }

        EXPECT_EQ(oxt_count_after, 0);
    }

    // Case 2: Internal residue (should not have OXT, function should be no-op)
    {
        // Target: chain A, residue 10
        gemmi::Residue* target_res = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10)
                            target_res = &res;

        ASSERT_NE(target_res, nullptr) << "residue //A/10 not found";

        int oxt_count_before = 0;
        for (const auto& atom : target_res->atoms) {
            std::string name = atom.name;
            while (!name.empty() && name.back() == ' ')
                name.pop_back();
            if (name == "OXT")
                oxt_count_before++;
        }

        EXPECT_EQ(oxt_count_before, 0);

        // Call the function
        coot::remove_OXT_internal_gemmi(target_res);

        int oxt_count_after = 0;
        for (const auto& atom : target_res->atoms) {
            std::string name = atom.name;
            while (!name.empty() && name.back() == ' ')
                name.pop_back();
            if (name == "OXT")
                oxt_count_after++;
        }

        EXPECT_EQ(oxt_count_after, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}