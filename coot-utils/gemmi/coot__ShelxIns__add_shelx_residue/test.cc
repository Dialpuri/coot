#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ShelxAddResidue) {
    coot::ShelxIns si;

    // Case 1: non-standard residue (should setHet=1)
    std::vector<gemmi::Atom> atom_vector;

    gemmi::Atom at1;
    at1.name = "C  ";
    at1.pos.x = 1.0;
    at1.pos.y = 2.0;
    at1.pos.z = 3.0;

    gemmi::Atom at2;
    at2.name = "N  ";
    at2.pos.x = 4.0;
    at2.pos.y = 5.0;
    at2.pos.z = 6.0;

    atom_vector.push_back(at1);
    atom_vector.push_back(at2);

    std::string res_name1 = "XYZ";
    int res_no1 = 10;

    gemmi::Residue* res1 = coot::add_shelx_residue_gemmi(si, atom_vector, res_name1, res_no1);

    EXPECT_TRUE(res1 != nullptr);
    EXPECT_EQ(res1->name, res_name1);
    EXPECT_EQ(res1->seqid.num.value, res_no1);
    EXPECT_EQ(atom_vector.size(), 2);

    // First atom should have occ=0 (Het=1 in original)
    EXPECT_EQ(atom_vector[0].occ, 0);

    // Case 2: standard residue (should NOT setHet=1)
    gemmi::Atom at3;
    at3.name = "CA ";
    at3.pos.x = 7.0;
    at3.pos.y = 8.0;
    at3.pos.z = 9.0;

    std::vector<gemmi::Atom> atom_vector2;
    atom_vector2.push_back(at3);

    std::string res_name2 = "ALA";
    int res_no2 = 42;

    gemmi::Residue* res2 = coot::add_shelx_residue_gemmi(si, atom_vector2, res_name2, res_no2);

    EXPECT_TRUE(res2 != nullptr);
    EXPECT_EQ(res2->name, res_name2);
    EXPECT_EQ(res2->seqid.num.value, res_no2);
    EXPECT_EQ(atom_vector2.size(), 1);

    // First atom should have occ=1 (Het=0 in original)
    EXPECT_EQ(atom_vector2[0].occ, 1);

    // Cleanup
    delete res1;
    delete res2;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}