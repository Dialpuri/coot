#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, add_energy_lib_bonds) {
    coot::energy_lib_t energy_lib;

    // Case 1: Valid bond — add one bond via a manually constructed cif::Loop
    {
        int bonds_before = static_cast<int>(energy_lib.bonds.size());

        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "type", "const", "length", "value_esd"});
        loop.add_row({"CA", "CB", "normal", "420.0", "1.54", "0.02"});

        EXPECT_EQ(bonds_before, 0);

        coot::add_energy_lib_bonds_gemmi(energy_lib, loop);

        int bonds_after = static_cast<int>(energy_lib.bonds.size());
        EXPECT_EQ(bonds_after, 1);
        EXPECT_EQ(bonds_after - bonds_before, 1);

        const coot::energy_lib_bond& b = energy_lib.bonds.back();
        EXPECT_EQ(b.atom_type_1, "CA");
        EXPECT_EQ(b.atom_type_2, "CB");
        EXPECT_EQ(b.type, "normal");
        EXPECT_NEAR(b.spring_constant, 420.0, 1e-4);
        EXPECT_NEAR(b.length, 1.54, 1e-4);
        EXPECT_NEAR(b.esd, 0.02, 1e-4);
    }

    // Case 2: delocal type gets converted to deloc
    {
        int bonds_before = static_cast<int>(energy_lib.bonds.size());

        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "type", "const", "length", "value_esd"});
        loop.add_row({"C", "O", "delocal", "500.0", "1.23", "0.01"});

        EXPECT_EQ(bonds_before, 1);

        coot::add_energy_lib_bonds_gemmi(energy_lib, loop);

        int bonds_after = static_cast<int>(energy_lib.bonds.size());
        EXPECT_EQ(bonds_after, 2);
        EXPECT_EQ(bonds_after - bonds_before, 1);

        const coot::energy_lib_bond& b = energy_lib.bonds.back();
        EXPECT_EQ(b.type, "deloc");
    }

    // Case 3: Missing const and value_esd → defaults applied (empty strings)
    {
        int bonds_before = static_cast<int>(energy_lib.bonds.size());

        gemmi::cif::Loop loop({"atom_type_1", "atom_type_2", "type", "const", "length", "value_esd"});
        loop.add_row({"N", "CA", "normal", "", "1.45", ""});

        EXPECT_EQ(bonds_before, 2);

        coot::add_energy_lib_bonds_gemmi(energy_lib, loop);

        int bonds_after = static_cast<int>(energy_lib.bonds.size());
        EXPECT_EQ(bonds_after, 3);
        EXPECT_EQ(bonds_after - bonds_before, 1);

        const coot::energy_lib_bond& b = energy_lib.bonds.back();
        EXPECT_NEAR(b.spring_constant, 420.0, 1e-4);
        EXPECT_NEAR(b.esd, 0.02, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}