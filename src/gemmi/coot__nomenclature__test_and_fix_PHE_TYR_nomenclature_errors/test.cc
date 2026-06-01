#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

TEST(OracleTest, test_and_fix_PHE_TYR_nomenclature_errors) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Case 1: Normal PHE residue with apply_swap=false (chi2 in range → no swap)
    {
        std::string cid = "//A/4";
        gemmi::Residue* res = coot::cid_to_residue_gemmi(cid, model);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(std::string(res->name), "PHE");

        // Find the chain pointer for this residue
        gemmi::Chain* ch = nullptr;
        for (auto& c : model.chains) {
            for (auto& r : c.residues) {
                if (&r == res) { ch = &c; break; }
            }
        }
        ASSERT_NE(ch, nullptr);

        gemmi::CRA cra{ch, res, nullptr};
        int result = coot::nomenclature::test_and_fix_PHE_TYR_nomenclature_errors_gemmi(cra, false);
        EXPECT_EQ(result, 1);
    }

    // Case 2: Normal TYR residue with apply_swap=true (chi2 in range → no swap)
    {
        std::string cid = "//A/15";
        gemmi::Residue* res = coot::cid_to_residue_gemmi(cid, model);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(std::string(res->name), "TYR");

        // Find the chain pointer for this residue
        gemmi::Chain* ch = nullptr;
        for (auto& c : model.chains) {
            for (auto& r : c.residues) {
                if (&r == res) { ch = &c; break; }
            }
        }
        ASSERT_NE(ch, nullptr);

        gemmi::CRA cra{ch, res, nullptr};
        int result = coot::nomenclature::test_and_fix_PHE_TYR_nomenclature_errors_gemmi(cra, true);
        EXPECT_EQ(result, 0);
    }

    // Case 3: Non-PHE/TYR residue (function finds no CD1 → returns 0)
    {
        std::string cid = "//A/1";
        gemmi::Residue* res = coot::cid_to_residue_gemmi(cid, model);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(std::string(res->name), "MET");

        // Find the chain pointer for this residue
        gemmi::Chain* ch = nullptr;
        for (auto& c : model.chains) {
            for (auto& r : c.residues) {
                if (&r == res) { ch = &c; break; }
            }
        }
        ASSERT_NE(ch, nullptr);

        gemmi::CRA cra{ch, res, nullptr};
        int result = coot::nomenclature::test_and_fix_PHE_TYR_nomenclature_errors_gemmi(cra, true);
        EXPECT_EQ(result, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}