#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

TEST(OracleTestGemmi, delete_residue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GT(st.models.size(), 0u);

    gemmi::Model& model = st.models[0];

    unsigned int atoms_before = static_cast<unsigned int>(gemmi::count_atom_sites(model));

    // Case 1: delete a real residue — chain A, residue 1
    {
        coot::residue_spec_t spec("A", 1, "");
        int result = coot::delete_residue_gemmi(st, model, spec);
        unsigned int atoms_after = static_cast<unsigned int>(gemmi::count_atom_sites(model));
        EXPECT_EQ(result, 1);
        EXPECT_EQ(atoms_before, 2107u);
        EXPECT_EQ(atoms_after, 2099u);
    }

    // Case 2: try to delete a non-existent residue (residue number too high)
    {
        coot::residue_spec_t spec("A", 9999, "");
        unsigned int atoms_before2 = static_cast<unsigned int>(gemmi::count_atom_sites(model));
        int result = coot::delete_residue_gemmi(st, model, spec);
        unsigned int atoms_after2 = static_cast<unsigned int>(gemmi::count_atom_sites(model));
        EXPECT_EQ(result, 0);
        EXPECT_EQ(atoms_before2, 2099u);
        EXPECT_EQ(atoms_after2, 2099u);
    }

    // Case 3: try to delete from a non-existent chain
    {
        coot::residue_spec_t spec("Z", 1, "");
        unsigned int atoms_before3 = static_cast<unsigned int>(gemmi::count_atom_sites(model));
        int result = coot::delete_residue_gemmi(st, model, spec);
        unsigned int atoms_after3 = static_cast<unsigned int>(gemmi::count_atom_sites(model));
        EXPECT_EQ(result, 0);
        EXPECT_EQ(atoms_before3, 2099u);
        EXPECT_EQ(atoms_after3, 2099u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}