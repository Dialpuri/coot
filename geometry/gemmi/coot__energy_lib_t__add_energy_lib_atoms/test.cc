#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, add_energy_lib_atoms) {
    coot::energy_lib_t energy_lib;

    // Case 1: Single atom (Oxygen acceptor)
    {
        std::size_t initial_size = energy_lib.atom_map.size();
        EXPECT_EQ(initial_size, 0u);

        gemmi::cif::Loop loop;
        loop.tags.push_back("type");
        loop.tags.push_back("weight");
        loop.tags.push_back("hb_type");
        loop.tags.push_back("vdw_radius");
        loop.tags.push_back("vdwh_radius");
        loop.tags.push_back("ion_radius");
        loop.tags.push_back("element");
        loop.tags.push_back("valency");
        loop.tags.push_back("sp");

        loop.add_row({ "O", "17.0", "A", "1.40", "1.20", "1.35", "O", "2", "2" });

        coot::add_energy_lib_atoms_gemmi(energy_lib, loop);

        std::size_t after_size = energy_lib.atom_map.size();
        EXPECT_EQ(after_size, 1u);

        EXPECT_TRUE(energy_lib.atom_map.count("O") > 0);
        coot::energy_lib_atom &at = energy_lib.atom_map["O"];
        EXPECT_EQ(at.type, "O");
        EXPECT_EQ(at.element, "O");
        EXPECT_NEAR(at.weight, 17.0, 1e-4);
        EXPECT_NEAR(at.vdw_radius, 1.4, 1e-4);
        EXPECT_EQ(static_cast<int>(at.hb_type), 2);
    }

    // Case 2: Two atoms in one loop (Nitrogen donor + Carbon neither)
    {
        gemmi::cif::Loop loop2;
        loop2.tags.push_back("type");
        loop2.tags.push_back("weight");
        loop2.tags.push_back("hb_type");
        loop2.tags.push_back("vdw_radius");
        loop2.tags.push_back("vdwh_radius");
        loop2.tags.push_back("ion_radius");
        loop2.tags.push_back("element");
        loop2.tags.push_back("valency");
        loop2.tags.push_back("sp");

        loop2.add_row({ "N", "14.0", "D", "1.50", "1.30", "1.25", "N", "3", "3" });
        loop2.add_row({ "C", "12.0", "N", "1.70", "1.50", "1.45", "C", "4", "3" });

        std::size_t before_size2 = energy_lib.atom_map.size();
        EXPECT_EQ(before_size2, 1u);

        coot::add_energy_lib_atoms_gemmi(energy_lib, loop2);

        std::size_t after_size2 = energy_lib.atom_map.size();
        EXPECT_EQ(after_size2, 3u);

        EXPECT_TRUE(energy_lib.atom_map.count("N") > 0);
        coot::energy_lib_atom &nat = energy_lib.atom_map["N"];
        EXPECT_EQ(nat.type, "N");
        EXPECT_EQ(static_cast<int>(nat.hb_type), 1);

        EXPECT_TRUE(energy_lib.atom_map.count("C") > 0);
        coot::energy_lib_atom &cat = energy_lib.atom_map["C"];
        EXPECT_EQ(cat.type, "C");
        EXPECT_EQ(static_cast<int>(cat.hb_type), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}