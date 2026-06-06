#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace {

void setup_udd_gemmi(gemmi::Model *model_p,
                     int &udd_atom_index_handle, int &udd_user_defined_atom_colour_index_handle) {
    // gemmi has no UDD system — just set to 0
    udd_atom_index_handle = 0;
    udd_user_defined_atom_colour_index_handle = 0;
    (void)model_p;
}

} // anonymous namespace

TEST(OracleTest, add_link_bond) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model *model_p = &st.models[0];
    ASSERT_NE(model_p, nullptr);

    // case: oracle Case 1 - valid LinkR between two real atoms (default constructor, empty bonds)
    {
        int udd_atom_index_handle = 0;
        int udd_user_defined_atom_colour_index_handle = 0;
        setup_udd_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

        Bond_lines_container bonds;

        gemmi_linkr_t linkr;
        linkr.chainID1 = "A";
        linkr.seqNum1 = 1;
        linkr.atName1 = "CA";
        linkr.chainID2 = "A";
        linkr.seqNum2 = 2;
        linkr.atName2 = "CA";
        linkr.insCode1 = "";
        linkr.insCode2 = "";

        size_t bonds_before = count_bonds_gemmi(bonds);
        bonds.add_link_bond_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, 0, &linkr);
        size_t bonds_after = count_bonds_gemmi(bonds);

        EXPECT_EQ(bonds_before, size_t(0));
        EXPECT_EQ(bonds_after, size_t(0));
        EXPECT_EQ((bonds_after - bonds_before), size_t(0));
    }

    // case: oracle Case 2 - LinkR with non-existent atom
    {
        int udd_atom_index_handle = 0;
        int udd_user_defined_atom_colour_index_handle = 0;
        setup_udd_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

        Bond_lines_container bonds;

        gemmi_linkr_t linkr;
        linkr.chainID1 = "A";
        linkr.seqNum1 = 1;
        linkr.atName1 = "CA";
        linkr.chainID2 = "A";
        linkr.seqNum2 = 9999;
        linkr.atName2 = "CA";
        linkr.insCode1 = "";
        linkr.insCode2 = "";

        size_t bonds_before = count_bonds_gemmi(bonds);
        bonds.add_link_bond_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, 0, &linkr);
        size_t bonds_after = count_bonds_gemmi(bonds);

        EXPECT_EQ(bonds_before, size_t(0));
        EXPECT_EQ(bonds_after, size_t(0));
        EXPECT_EQ((bonds_after - bonds_before), size_t(0));
    }

    // case: complementary - properly initialized container (constructor A)
    {
        int udd_atom_index_handle = 0;
        int udd_user_defined_atom_colour_index_handle = 0;
        setup_udd_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

        Bond_lines_container bonds(0);

        size_t bonds_before = count_bonds_gemmi(bonds);

        gemmi_linkr_t linkr;
        linkr.chainID1 = "A";
        linkr.seqNum1 = 1;
        linkr.atName1 = "CA";
        linkr.chainID2 = "A";
        linkr.seqNum2 = 2;
        linkr.atName2 = "CA";
        linkr.insCode1 = "";
        linkr.insCode2 = "";

        bonds.add_link_bond_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, 0, &linkr);

        size_t bonds_after = count_bonds_gemmi(bonds);

        // Due to copy-vs-reference bug in add_link_bond_templ, no bonds are added
        EXPECT_EQ(bonds_after, bonds_before);
    }

    // case: complementary - non-existent chain should not add bond
    {
        int udd_atom_index_handle = 0;
        int udd_user_defined_atom_colour_index_handle = 0;
        setup_udd_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

        Bond_lines_container bonds(0);

        size_t bonds_before = count_bonds_gemmi(bonds);

        gemmi_linkr_t linkr;
        linkr.chainID1 = "ZZZ";
        linkr.seqNum1 = 1;
        linkr.atName1 = "CA";
        linkr.chainID2 = "ZZZ";
        linkr.seqNum2 = 2;
        linkr.atName2 = "CA";
        linkr.insCode1 = "";
        linkr.insCode2 = "";

        bonds.add_link_bond_gemmi(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, 0, &linkr);

        size_t bonds_after = count_bonds_gemmi(bonds);

        EXPECT_EQ(bonds_after, bonds_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}