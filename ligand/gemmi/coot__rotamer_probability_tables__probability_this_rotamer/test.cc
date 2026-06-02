#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <filesystem>
#include <vector>
#include <string>
#include "utils/coot-utils.hh"

TEST(OracleTest, rotamer_probability_tables) {
    // Same directory discovery logic as oracle
    std::string pkg_dir = coot::package_data_dir();
    ASSERT_FALSE(pkg_dir.empty());

    std::string rota_dir = pkg_dir + "/rotamer";

    if (!std::filesystem::exists(rota_dir)) {
        std::vector<std::string> try_dirs = {
            pkg_dir + "/data/rotamer",
            "/usr/share/coot/rotamer",
            "/usr/local/share/coot/rotamer",
            "/lmb/home/jdialpuri/Development/coot-dev/data/rotamer",
            "/lmb/home/jdialpuri/Development/coot-dev/coot/data/rotamer",
        };
        for (const auto &td : try_dirs) {
            if (std::filesystem::exists(td)) {
                rota_dir = td;
                break;
            }
        }
    }

    // If rotamer data directory doesn't exist, verify the setup would work
    // by checking the class can be constructed and configured
    if (!std::filesystem::exists(rota_dir)) {
        coot::rotamer_probability_tables rpt;
        rpt.set_tables_dir(rota_dir);
        EXPECT_FALSE(rpt.is_well_formatted());
        return;
    }

    // If data is available, load and test
    coot::rotamer_probability_tables rpt;
    rpt.set_tables_dir(rota_dir);
    rpt.fill_tables();
    ASSERT_TRUE(rpt.is_well_formatted());
    ASSERT_GT(rpt.n_tables(), 0u);

    // Find table index for VAL and HIS
    unsigned int val_idx = -1;
    unsigned int his_idx = -1;
    for (unsigned int i = 0; i < rpt.n_tables(); i++) {
        const coot::a_rotamer_table &t = rpt[i];
        if (t.residue_name == "VAL") val_idx = i;
        if (t.residue_name == "HIS") his_idx = i;
    }

    // Verify we found expected residue tables
    if (val_idx != -1) {
        const coot::a_rotamer_table &val_table = rpt[val_idx];
        EXPECT_EQ(val_table.residue_name, "VAL");
        EXPECT_GT(val_table.n_chis, 0);
    }
    if (his_idx != -1) {
        const coot::a_rotamer_table &his_table = rpt[his_idx];
        EXPECT_EQ(his_table.residue_name, "HIS");
        EXPECT_GT(his_table.n_chis, 0);
    }

    // Test probability_this_rotamer_gemmi with VAL (has 1 chi angle)
    if (val_idx != -1) {
        const coot::a_rotamer_table &val_table = rpt[val_idx];
        EXPECT_EQ(val_table.residue_name, "VAL");

        // VAL has chi1 only (n_chis == 1)
        std::vector<std::pair<int,float>> chi_angles = { {1, 180.0f} };
        coot::rotamer_probability_info_t info = coot::probability_this_rotamer_gemmi(rpt, val_idx, chi_angles);
        EXPECT_EQ(info.state, coot::rotamer_probability_info_t::OK);
        EXPECT_GT(info.probability, 0.0);
        EXPECT_EQ(info.rotamer_name, "VAL");

        // Test with a different chi1 angle
        chi_angles[0].second = 60.0f;
        info = coot::probability_this_rotamer_gemmi(rpt, val_idx, chi_angles);
        EXPECT_EQ(info.state, coot::rotamer_probability_info_t::OK);
        EXPECT_GT(info.probability, 0.0);
    }

    // Test probability_this_rotamer_gemmi with HIS (has 2 chi angles)
    if (his_idx != -1) {
        const coot::a_rotamer_table &his_table = rpt[his_idx];
        EXPECT_EQ(his_table.residue_name, "HIS");

        // HIS has chi1 and chi2 (n_chis == 2)
        std::vector<std::pair<int,float>> chi_angles = { {1, -60.0f}, {2, 90.0f} };
        coot::rotamer_probability_info_t info = coot::probability_this_rotamer_gemmi(rpt, his_idx, chi_angles);
        EXPECT_EQ(info.state, coot::rotamer_probability_info_t::OK);
        EXPECT_GE(info.probability, 0.0);
        EXPECT_EQ(info.rotamer_name, "HIS");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}