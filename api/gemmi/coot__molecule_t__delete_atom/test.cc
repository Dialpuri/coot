#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeleteAtom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    unsigned int atoms_before = gemmi::count_atom_sites(st);
    EXPECT_EQ(atoms_before, 2107u);

    {
        coot::atom_spec_t atom_spec;
        atom_spec.chain_id = "A";
        atom_spec.res_no = 12;
        atom_spec.ins_code = "";
        atom_spec.atom_name = " OE1 ";
        atom_spec.alt_conf = "";

        std::pair<int, unsigned int> result = coot::delete_atom_gemmi(st, atom_spec);
        int was_deleted = result.first;
        unsigned int atoms_after = result.second;

        EXPECT_EQ(was_deleted, 0);
        EXPECT_EQ(atoms_after, 2107u);
    }

    {
        coot::atom_spec_t atom_spec;
        atom_spec.chain_id = "A";
        atom_spec.res_no = 14;
        atom_spec.ins_code = "";
        atom_spec.atom_name = " OG1 ";
        atom_spec.alt_conf = "";

        std::pair<int, unsigned int> result = coot::delete_atom_gemmi(st, atom_spec);
        int was_deleted = result.first;
        unsigned int atoms_after = result.second;

        EXPECT_EQ(was_deleted, 0);
        EXPECT_EQ(atoms_after, 2107u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}