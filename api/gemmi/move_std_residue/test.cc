#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_standard_residue_instance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__deep_copy_this_residue_old_style/gemmi/function.hh"

TEST(OracleTest, move_std_residue) {
    // Case 1: Normal call with valid reference residue
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        auto& model = st.models[0];

        auto& chain = model.chains[0];

        ASSERT_GT(chain.residues.size(), 0);

        auto& reference_residue = chain.residues[0];

        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(reference_residue.seqid.num.value, 1);
        EXPECT_EQ(reference_residue.name, "MET");

        // Get a standard residue instance (ALA)
        gemmi::Residue* std_res = coot::util::get_standard_residue_instance_gemmi("ALA");
        ASSERT_NE(std_res, nullptr);

        // Build CRA for the standard residue (no parent chain)
        gemmi::CRA std_cra{nullptr, std_res, nullptr};

        // Deep copy it so we can move it
        gemmi::Residue* std_res_copy = coot::deep_copy_this_residue_old_style_gemmi(std_cra, "", 1, false);
        ASSERT_NE(std_res_copy, nullptr);
        ASSERT_FALSE(std_res_copy->atoms.empty());

        // Print some atoms from the moving residue before the call
        EXPECT_EQ(std_res_copy->atoms.size(), 5);

        // First atom before the call
        const gemmi::Atom& first_atom = std_res_copy->atoms[0];
        EXPECT_EQ(first_atom.name, "N");
        EXPECT_NEAR(first_atom.pos.x, 0.824, 1e-4);
        EXPECT_NEAR(first_atom.pos.y, 0.0, 1e-4);
        EXPECT_NEAR(first_atom.pos.z, 1.203, 1e-4);

        // Build CRA for moving residue (no parent chain for a standard residue)
        gemmi::CRA cra{nullptr, std_res_copy, nullptr};

        // Call move_std_residue
        int result = move_std_residue_gemmi(cra, &reference_residue);

        EXPECT_EQ(result, 1);

        // Print atoms after the call to verify coordinates changed
        ASSERT_GT(std_res_copy->atoms.size(), 0);
        const gemmi::Atom& first_atom_after = std_res_copy->atoms[0];
        EXPECT_EQ(first_atom_after.name, "N");
        EXPECT_NEAR(first_atom_after.pos.x, 12.3529, 1e-4);
        EXPECT_NEAR(first_atom_after.pos.y, 2.79091, 1e-4);
        EXPECT_NEAR(first_atom_after.pos.z, 74.9254, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}