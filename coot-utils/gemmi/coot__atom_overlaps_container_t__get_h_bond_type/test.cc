#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_h_bond_type_hardcoded) {
    coot::protein_geometry geom;
    geom.init_standard();

    coot::atom_overlaps_container_t container(&geom);

    // Create a gemmi Residue with name "MET"
    gemmi::Residue res;
    res.name = "MET";

    // Create a gemmi Atom with 4-char padded name " N  " (as in PDB format)
    gemmi::Atom atom;
    atom.name = " N  ";

    {
        coot::hb_t type = container.get_h_bond_type_gemmi(atom, res);
        EXPECT_EQ(type, coot::hb_t::HB_DONOR);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}