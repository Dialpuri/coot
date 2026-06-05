#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, add_energy_lib_torsions) {
    coot::energy_lib_t energy_lib;

    // Case 1: Valid torsion data — manually constructed Loop with 2 rows
    {
        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "atom_type_3", "atom_type_4",
                               "const", "angle", "period"});
        loop.add_row({"C", "C", "C", "C", "2.5", "180.0", "2"});
        loop.add_row({"CA", "CA", "C", "O", "0.8", "0.0", "2"});

        EXPECT_EQ(loop.length(), 2u);
        EXPECT_EQ(energy_lib.torsions.size(), 0u);

        coot::add_energy_lib_torsions_gemmi(energy_lib, loop);

        EXPECT_EQ(energy_lib.torsions.size(), 2u);

        coot::energy_lib_torsion &t0 = energy_lib.torsions[0];
        EXPECT_EQ(t0.atom_type_1, "C");
        EXPECT_EQ(t0.atom_type_2, "C");
        EXPECT_EQ(t0.atom_type_3, "C");
        EXPECT_EQ(t0.atom_type_4, "C");
        EXPECT_NEAR(t0.spring_constant, 2.5, 1e-4);
        EXPECT_NEAR(t0.angle, 180.0, 1e-4);
        EXPECT_EQ(t0.period, 2);

        coot::energy_lib_torsion &t1 = energy_lib.torsions[1];
        EXPECT_EQ(t1.atom_type_1, "CA");
        EXPECT_EQ(t1.atom_type_2, "CA");
        EXPECT_EQ(t1.atom_type_3, "C");
        EXPECT_EQ(t1.atom_type_4, "O");
        EXPECT_NEAR(t1.spring_constant, 0.8, 1e-4);
        EXPECT_NEAR(t1.angle, 0.0, 1e-4);
        EXPECT_EQ(t1.period, 2);
    }

    // Case 2: Empty loop — should add nothing
    {
        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "atom_type_3", "atom_type_4",
                               "const", "angle", "period"});
        // No rows added

        EXPECT_EQ(loop.length(), 0u);
        EXPECT_EQ(energy_lib.torsions.size(), 2u);

        coot::add_energy_lib_torsions_gemmi(energy_lib, loop);

        EXPECT_EQ(energy_lib.torsions.size(), 2u);
    }

    // Case 3: Row with bad period value ("!!" instead of integer) — should be skipped by guard
    {
        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "atom_type_3", "atom_type_4",
                               "const", "angle", "period"});
        loop.add_row({"N", "CA", "C", "O", "1.0", "180.0", "!!"});

        EXPECT_EQ(loop.length(), 1u);
        EXPECT_EQ(energy_lib.torsions.size(), 2u);

        coot::add_energy_lib_torsions_gemmi(energy_lib, loop);

        EXPECT_EQ(energy_lib.torsions.size(), 2u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}