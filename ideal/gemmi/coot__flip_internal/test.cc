#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, flip_internal) {
    // --- setup ---
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find residues A:1 and A:2
    gemmi::Residue *res1 = nullptr, *res2 = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != "A") continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1 && res.seqid.icode == ' ') res1 = &res;
                if (res.seqid.num.value == 2 && res.seqid.icode == ' ') res2 = &res;
            }
        }
    }
    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    // Find atoms by name — gemmi stores names without padding
    gemmi::Atom *ca1 = nullptr, *ca2 = nullptr, *c1 = nullptr, *o1 = nullptr, *n_atm = nullptr;
    for (auto& atom : res1->atoms) {
        std::string an = atom.name;
        if (an == "CA" || an == " CA ") ca1 = &atom;
        if (an == "C"  || an == " C  ") c1 = &atom;
        if (an == "O"  || an == " O  ") o1 = &atom;
    }
    for (auto& atom : res2->atoms) {
        std::string an = atom.name;
        if (an == "CA" || an == " CA ") ca2 = &atom;
        if (an == "N"  || an == " N  ") n_atm = &atom;
    }
    ASSERT_NE(ca1, nullptr);
    ASSERT_NE(ca2, nullptr);
    ASSERT_NE(c1, nullptr);
    ASSERT_NE(o1, nullptr);
    ASSERT_NE(n_atm, nullptr);

    // Build CAS vector once, used by both cases
    std::vector<gemmi::Vec3> cas(2);
    cas[0] = gemmi::Vec3(ca1->pos.x, ca1->pos.y, ca1->pos.z);
    cas[1] = gemmi::Vec3(ca2->pos.x, ca2->pos.y, ca2->pos.z);

    // Case 1: normal flip
    {
        std::cout << "INPUT  ca1_orig: " << ca1->pos.x << " " << ca1->pos.y << " " << ca1->pos.z << std::endl;
        std::cout << "INPUT  ca2_orig: " << ca2->pos.x << " " << ca2->pos.y << " " << ca2->pos.z << std::endl;
        std::cout << "INPUT  c1_orig: " << c1->pos.x << " " << c1->pos.y << " " << c1->pos.z << std::endl;
        std::cout << "INPUT  o1_orig: " << o1->pos.x << " " << o1->pos.y << " " << o1->pos.z << std::endl;
        std::cout << "INPUT  n2_orig: " << n_atm->pos.x << " " << n_atm->pos.y << " " << n_atm->pos.z << std::endl;

        std::vector<gemmi::Vec3> flipping_atoms;
        flipping_atoms.push_back(c1->pos);
        flipping_atoms.push_back(o1->pos);
        flipping_atoms.push_back(n_atm->pos);

        std::vector<gemmi::Vec3> flipped = coot::flip_internal_gemmi(cas, flipping_atoms);

        // OUTPUT flipped_size: 3
        EXPECT_EQ(flipped.size(), 3);

        // OUTPUT c1_flipped: 12.9372 3.96125 76.1454
        EXPECT_NEAR(flipped[0].x, 12.9372, 1e-4);
        EXPECT_NEAR(flipped[0].y, 3.96125, 1e-4);
        EXPECT_NEAR(flipped[0].z, 76.1454, 1e-4);

        // OUTPUT o1_flipped: 11.7682 3.72864 76.4452
        EXPECT_NEAR(flipped[1].x, 11.7682, 1e-4);
        EXPECT_NEAR(flipped[1].y, 3.72864, 1e-4);
        EXPECT_NEAR(flipped[1].z, 76.4452, 1e-4);

        // OUTPUT n2_flipped: 13.8592 4.37744 77.0262
        EXPECT_NEAR(flipped[2].x, 13.8592, 1e-4);
        EXPECT_NEAR(flipped[2].y, 4.37744, 1e-4);
        EXPECT_NEAR(flipped[2].z, 77.0262, 1e-4);

        // OUTPUT c1_translation_magnitude: 1.08467
        double c1_dx = flipped[0].x - c1->pos.x;
        double c1_dy = flipped[0].y - c1->pos.y;
        double c1_dz = flipped[0].z - c1->pos.z;
        double dist = std::sqrt(c1_dx*c1_dx + c1_dy*c1_dy + c1_dz*c1_dz);
        EXPECT_NEAR(dist, 1.08467, 1e-4);
    }

    // Case 2: empty atom vector
    {
        std::vector<gemmi::Vec3> empty_flipped = coot::flip_internal_gemmi(cas, std::vector<gemmi::Vec3>());
        // OUTPUT empty_atoms_size: 0
        EXPECT_EQ(empty_flipped.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}