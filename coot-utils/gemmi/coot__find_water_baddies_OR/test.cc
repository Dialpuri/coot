#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <vector>
#include <string>
#include <unistd.h>

#include "api/molecules-container.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, find_water_baddies_OR) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_GT(st.models.size(), 0u);

    // Load MTZ for map
    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);
    EXPECT_GE(imol_map, 0);

    clipper::Xmap<float> xmap = mc.get_xmap(imol_map);
    float map_sigma = mc.get_map_rmsd_approx(imol_map);

    // Parameters
    float b_factor_lim = 80.0f;
    float outlier_sigma_level = 1.5f;
    float min_dist = 1.5f;
    float max_dist = 4.0f;
    short int ignore_part_occ_contact_flag = 0;
    short int ignore_zero_occ_flag = 0;

    // Collect all non-hydrogen atoms as selected_atoms
    std::vector<gemmi::Atom*> selected_atoms;
    for (gemmi::Chain& chain : st.models[0].chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                if (!atom.element.is_hydrogen()) {
                    selected_atoms.push_back(&atom);
                }
            }
        }
    }

    // Case 1: default parameters — no baddies found
    {
        std::vector<coot::atom_spec_t> result = coot::find_water_baddies_OR_gemmi(
            st.models[0],
            1,
            selected_atoms,
            b_factor_lim,
            xmap,
            map_sigma,
            outlier_sigma_level,
            min_dist, max_dist,
            ignore_part_occ_contact_flag,
            ignore_zero_occ_flag
        );

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: Very low b_factor_lim to flag all waters with tempFactor > threshold
    {
        std::vector<coot::atom_spec_t> result2 = coot::find_water_baddies_OR_gemmi(
            st.models[0],
            1,
            selected_atoms,
            0.5f,
            xmap,
            map_sigma,
            -100.0f,  // disable map sigma test
            min_dist, max_dist,
            ignore_part_occ_contact_flag,
            ignore_zero_occ_flag
        );

        EXPECT_EQ(result2.size(), 1u);
        std::cout << "DEBUG format: [" << result2[0].format() << "]" << std::endl;
        EXPECT_EQ(result2[0].format(), std::string("[spec: model 1 \"B\"    1 \"\" \" O  \" \"\"]"));
    }

    // Case 3: disable distance tests via negative values
    {
        std::vector<coot::atom_spec_t> result3 = coot::find_water_baddies_OR_gemmi(
            st.models[0],
            1,
            selected_atoms,
            0.5f,
            xmap,
            map_sigma,
            -100.0f,  // disable map sigma test
            -1.0f, -1.0f,  // disable distance tests
            ignore_part_occ_contact_flag,
            ignore_zero_occ_flag
        );

        EXPECT_EQ(result3.size(), 1u);
        EXPECT_EQ(result3[0].format(), std::string("[spec: model 1 \"B\"    1 \"\" \" O  \" \"\"]"));
    }
}
