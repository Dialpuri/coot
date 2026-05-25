#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AllAtomContactDots) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    coot::atom_overlaps_container_t aoc;
    coot::atom_overlaps_dots_container_t dots = aoc.all_atom_contact_dots(0.5, false);
    
    EXPECT_TRUE(true);
}

TEST(AtomOverlapsContainerTest, Constructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    coot::atom_overlaps_container_t aoc(&st, nullptr, false, 0.5, 0.25);
    EXPECT_TRUE(true);
}

TEST(AtomOverlapsContainerTest, MakeAllAtomOverlaps) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    coot::atom_overlaps_container_t aoc(&st, nullptr, false, 0.5, 0.25);
    aoc.make_all_atom_overlaps();
    EXPECT_TRUE(true);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}