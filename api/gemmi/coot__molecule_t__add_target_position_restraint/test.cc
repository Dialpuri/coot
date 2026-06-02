#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_target_position_restraint_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB load failed";

    std::vector<std::pair<gemmi::Atom *, clipper::Coord_orth>> atoms_with_position_restraints;

    // case: valid atom — add a new position restraint
    {
        std::string cid = "//A/10/CA";
        float pos_x = 10.0f, pos_y = 20.0f, pos_z = 30.0f;

        size_t before = atoms_with_position_restraints.size();

        coot::molecule_t::add_target_position_restraint_gemmi(cid, pos_x, pos_y, pos_z, st, atoms_with_position_restraints);

        size_t after = atoms_with_position_restraints.size();

        EXPECT_EQ(before, size_t(0));
        EXPECT_EQ(after, size_t(1));

        // Print the stored position to confirm it was set
        auto &restr = atoms_with_position_restraints.back();
        EXPECT_NEAR(restr.second.x(), 10.0f, 1e-4);
        EXPECT_NEAR(restr.second.y(), 20.0f, 1e-4);
        EXPECT_NEAR(restr.second.z(), 30.0f, 1e-4);
    }

    // case: same atom again — should update existing restraint (size stays the same)
    {
        std::string cid = "//A/10/CA";
        float pos_x = 50.0f, pos_y = 60.0f, pos_z = 70.0f;

        size_t before = atoms_with_position_restraints.size();

        coot::molecule_t::add_target_position_restraint_gemmi(cid, pos_x, pos_y, pos_z, st, atoms_with_position_restraints);

        size_t after = atoms_with_position_restraints.size();

        EXPECT_EQ(before, size_t(1));
        EXPECT_EQ(after, size_t(1));

        // Print the stored position — should be updated
        auto &restr = atoms_with_position_restraints.back();
        EXPECT_NEAR(restr.second.x(), 50.0f, 1e-4);
        EXPECT_NEAR(restr.second.y(), 60.0f, 1e-4);
        EXPECT_NEAR(restr.second.z(), 70.0f, 1e-4);
    }

    // case: invalid CID — nothing should change
    {
        std::string cid = "//A/9999/N";
        float pos_x = 100.0f, pos_y = 200.0f, pos_z = 300.0f;

        size_t before = atoms_with_position_restraints.size();

        coot::molecule_t::add_target_position_restraint_gemmi(cid, pos_x, pos_y, pos_z, st, atoms_with_position_restraints);

        size_t after = atoms_with_position_restraints.size();

        EXPECT_EQ(before, size_t(1));
        EXPECT_EQ(after, size_t(1));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}