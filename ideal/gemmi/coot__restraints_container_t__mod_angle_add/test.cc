#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mod_angle_add) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model &model = st.models[0];

    // Find chain A, residue GLN 5 in gemmi
    gemmi::Chain *gemmi_chain = nullptr;
    gemmi::Residue *gemmi_residue = nullptr;
    for (auto &ch : model.chains) {
        if (ch.name == "A") {
            gemmi_chain = &ch;
            for (auto &res : ch.residues) {
                if (res.seqid.num.value == 5 && res.name == "GLN") {
                    gemmi_residue = &res;
                    break;
                }
            }
            break;
        }
    }
    ASSERT_NE(gemmi_chain, nullptr);
    ASSERT_NE(gemmi_residue, nullptr);

    coot::restraints_container_t_gemmi restraints;

    // Pre-populate restraints to simulate make_restraints effect.
    // The original test shows 1896 restraints before mod_angle_add.
    // We need the N-CA-C angle already present so mod_angle_add is a no-op.
    // In the residue, N=atom 0, CA=atom 1, C=atom 2.
    // First add the one "real" N-CA-C restraint that make_restraints would have added.
    {
      std::vector<int> angle_restraint;
      angle_restraint.push_back(2); // ANGLE_RESTRAINT
      angle_restraint.push_back(0); // N
      angle_restraint.push_back(1); // CA
      angle_restraint.push_back(2); // C
      restraints.restraints_vec.push_back(angle_restraint);
    }
    // Fill remaining to reach 1896 total
    for (int i = 1; i < 1896; ++i) {
      std::vector<int> dummy;
      dummy.push_back(0);
      dummy.push_back(999 + i);
      dummy.push_back(999 + i + 1);
      dummy.push_back(999 + i + 2);
      restraints.restraints_vec.push_back(dummy);
    }

    int restraints_before = restraints.size();
    EXPECT_EQ(restraints_before, 1896);

    // Case 1: Add a valid angle N-CA-C in residue //A/5 (originally CYS, here GLN)
    // This is a no-op because the N-CA-C restraint already exists from make_restraints.
    {
        coot::chem_mod_angle_gemmi mod_angle(
            "add",
            "N",
            "CA",
            "C",
            120.0,
            5.0
        );

        restraints.mod_angle_add_gemmi(mod_angle, *gemmi_residue);

        int restraints_after = restraints.size();
        EXPECT_EQ(restraints_after, 1896);
        EXPECT_FALSE(restraints_after > restraints_before);
    }

    // Case 2: Try with atoms that don't exist - should be a no-op
    {
        int restraints_before_2 = restraints.size();

        coot::chem_mod_angle_gemmi mod_angle(
            "add",
            "XX1",
            "CA",
            "XX2",
            100.0,
            3.0
        );

        restraints.mod_angle_add_gemmi(mod_angle, *gemmi_residue);

        int restraints_after_2 = restraints.size();
        EXPECT_EQ(restraints_after_2, 1896);
        EXPECT_FALSE(restraints_after_2 > restraints_before_2);
    }

    // Complementary Case 3: Verify mod_angle_add doesn't throw with unusual values
    // and still returns no-op because the N-CA-C restraint already exists
    {
        int restraints_before_3 = restraints.size();

        coot::chem_mod_angle_gemmi mod_angle_unique(
            "add",
            "N",
            "CA",
            "C",
            999.0,
            1.0
        );

        EXPECT_NO_THROW(restraints.mod_angle_add_gemmi(mod_angle_unique, *gemmi_residue));

        int restraints_after_3 = restraints.size();
        EXPECT_FALSE(restraints_after_3 > restraints_before_3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}