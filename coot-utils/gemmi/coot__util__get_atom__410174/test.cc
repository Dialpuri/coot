#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(GetAtomTest, BasicFunctionality) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid atom - CA of residue 10 in chain A
    {
        gemmi::Atom *atom = coot::util::get_atom_gemmi("A", 10, "", "CA", "", st);
        if (atom) {
            std::cout << "atom_name: " << atom->name << std::endl;
            std::cout << "atom_coords: " << atom->pos.x << ", " << atom->pos.y << ", " << atom->pos.z << std::endl;
        }
        EXPECT_TRUE(atom != nullptr);
    }

    // case 2: valid atom - N of residue 1 in chain A
    {
        gemmi::Atom *atom = coot::util::get_atom_gemmi("A", 1, "", "N", "", st);
        if (atom) {
            std::cout << "atom_name: " << atom->name << std::endl;
        }
        EXPECT_TRUE(atom != nullptr);
    }

    // case 3: invalid atom - non-existent atom name
    {
        gemmi::Atom *atom = coot::util::get_atom_gemmi("A", 10, "", "XYZ", "", st);
        EXPECT_TRUE(atom == nullptr);
    }

    // case 4: invalid atom - non-existent residue
    {
        gemmi::Atom *atom = coot::util::get_atom_gemmi("A", 9999, "", "CA", "", st);
        EXPECT_TRUE(atom == nullptr);
    }

    // case 5: valid atom with altconf - O of residue 1 in chain A
    {
        gemmi::Atom *atom = coot::util::get_atom_gemmi("A", 1, "", "O", "", st);
        if (atom) {
            std::cout << "atom_name: " << atom->name << std::endl;
        }
        EXPECT_TRUE(atom != nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}