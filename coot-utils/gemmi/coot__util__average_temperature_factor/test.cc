#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, average_temperature_factor) {
    // Setup: load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get residue A/10
    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = model.find_chain("A");
    ASSERT_NE(chain, nullptr);

    const gemmi::Residue* res = nullptr;
    for (auto& r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    // Get atoms from the residue into a vector of pointers
    std::vector<const gemmi::Atom*> residue_atoms;
    for (auto& atom : res->atoms) {
        residue_atoms.push_back(&atom);
    }
    int n_residue_atoms = static_cast<int>(residue_atoms.size());

    EXPECT_EQ(n_residue_atoms, 8);

    // Case 1: No cutoffs — average of all atoms' tempFactors / n_atoms
    {
        float low_cutoff = 0.0f;
        float high_cutoff = 100.0f;
        short int apply_low = 0;
        short int apply_high = 0;

        float avg = coot::util::average_temperature_factor_gemmi(
            residue_atoms, n_residue_atoms,
            low_cutoff, high_cutoff,
            apply_low, apply_high);

        EXPECT_NEAR(avg, 40.2838f, 1e-4);
    }

    // Case 2: Apply low cutoff — only atoms with tempFactor > 30.0 are included
    {
        float low_cutoff = 30.0f;
        float high_cutoff = 0.0f;
        short int apply_low = 1;
        short int apply_high = 0;

        float avg = coot::util::average_temperature_factor_gemmi(
            residue_atoms, n_residue_atoms,
            low_cutoff, high_cutoff,
            apply_low, apply_high);

        EXPECT_NEAR(avg, 36.76f, 1e-4);

        // Verify individual tempFactors match oracle output
        EXPECT_NEAR(residue_atoms[0]->b_iso, 38.99f, 1e-3);
        EXPECT_NEAR(residue_atoms[1]->b_iso, 41.55f, 1e-3);
        EXPECT_NEAR(residue_atoms[2]->b_iso, 49.47f, 1e-3);
        EXPECT_NEAR(residue_atoms[3]->b_iso, 52.82f, 1e-3);
        EXPECT_NEAR(residue_atoms[4]->b_iso, 38.34f, 1e-3);
        EXPECT_NEAR(residue_atoms[5]->b_iso, 35.8f, 1e-3);
        EXPECT_NEAR(residue_atoms[6]->b_iso, 37.11f, 1e-3);
        EXPECT_NEAR(residue_atoms[7]->b_iso, 28.19f, 1e-3);
    }

    // Case 3: Edge case — zero atoms (should return 0.0)
    {
        std::vector<const gemmi::Atom*> empty_atoms;
        int n_empty = 0;
        float avg = coot::util::average_temperature_factor_gemmi(
            empty_atoms, n_empty,
            0.0f, 100.0f,
            0, 0);

        EXPECT_NEAR(avg, 0.0f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
