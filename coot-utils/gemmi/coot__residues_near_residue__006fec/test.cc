#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

inline std::string norm_ic(const std::string& ic) {
    return ic.empty() ? std::string(" ") : ic;
}

inline bool contains_residue(const std::vector<coot::residue_spec_t>& vec,
                             const coot::residue_spec_t& spec) {
    for (const auto& r : vec) {
        if (norm_ic(r.chain_id) == norm_ic(spec.chain_id) &&
            r.res_no == spec.res_no &&
            norm_ic(r.ins_code) == norm_ic(spec.ins_code)) {
            return true;
        }
    }
    return false;
}

TEST(OracleTest, ResiduesNearResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    {
        coot::residue_spec_t spec("A", 10, " ");
        float radius = 4.0f;
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(spec, st, radius);
        
        EXPECT_EQ(result.size(), 6);
        
        std::vector<coot::residue_spec_t> expected_specs = {
            coot::residue_spec_t("A", 18, " "),
            coot::residue_spec_t("A", 8, " "),
            coot::residue_spec_t("A", 11, " "),
            coot::residue_spec_t("A", 9, " "),
            coot::residue_spec_t("A", 19, " "),
            coot::residue_spec_t("A", 20, " ")
        };
        
        for (const auto& exp_spec : expected_specs) {
            EXPECT_TRUE(contains_residue(result, exp_spec))
                << "Expected residue not found: " << exp_spec.chain_id 
                << "/" << exp_spec.res_no << "/" << exp_spec.ins_code;
        }
    }

    {
        coot::residue_spec_t spec("A", 50, " ");
        float radius = 5.0f;
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(spec, st, radius);
        
        EXPECT_EQ(result.size(), 12);
    }

    {
        coot::residue_spec_t spec("A", 9999, " ");
        float radius = 4.0f;
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(spec, st, radius);
        
        EXPECT_EQ(result.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}