#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, add_o2_prime) {
    // Test case 1: Standard RNA residue with atoms at "expected" coordinates
    {
        gemmi::Residue res_copy;
        res_copy.name = "G";
        res_copy.seqid.num = 1;  // use int assignment, no insertion code
        // res_copy.seqid.icode = ' '; // default is space

        // Add C1', C2', C3' atoms with typical RNA coordinates
        gemmi::Atom atom_c1p;
        atom_c1p.name = " C1'";
        atom_c1p.element = gemmi::Element("C");
        atom_c1p.pos = gemmi::Position(6.853, -5.219, 1.725);
        atom_c1p.occ = 1.0f;
        atom_c1p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c1p);

        gemmi::Atom atom_c2p;
        atom_c2p.name = " C2'";
        atom_c2p.element = gemmi::Element("C");
        atom_c2p.pos = gemmi::Position(7.473, -4.888, 3.087);
        atom_c2p.occ = 1.0f;
        atom_c2p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c2p);

        gemmi::Atom atom_c3p;
        atom_c3p.name = " C3'";
        atom_c3p.element = gemmi::Element("C");
        atom_c3p.pos = gemmi::Position(6.678, -5.807, 4.033);
        atom_c3p.occ = 1.0f;
        atom_c3p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c3p);

        EXPECT_EQ(res_copy.atoms.size(), 3U);

        coot::ideal_rna::add_o2_prime_gemmi(&res_copy);

        EXPECT_EQ(res_copy.atoms.size(), 4U);

        bool has_o2p = false;
        for (const auto& atom : res_copy.atoms) {
            if (atom.name.find("O2'") != std::string::npos ||
                atom.name.find("O2*") != std::string::npos) {
                has_o2p = true;
                // Expected: (8.87, -5.158, 3.018) with some tolerance
                EXPECT_NEAR(atom.pos.x, 8.87, 0.01);
                EXPECT_NEAR(atom.pos.y, -5.158, 0.01);
                EXPECT_NEAR(atom.pos.z, 3.018, 0.01);
            }
        }
        EXPECT_TRUE(has_o2p);
    }

    // Test case 2: RNA residue with atoms at different coordinates (shifted)
    {
        gemmi::Residue res_copy;
        res_copy.name = "G";
        res_copy.seqid.num = 2;
        // res_copy.seqid.icode = ' ';

        // Add C1', C2', C3' atoms shifted by (10, 0, 0)
        gemmi::Atom atom_c1p;
        atom_c1p.name = " C1'";
        atom_c1p.element = gemmi::Element("C");
        atom_c1p.pos = gemmi::Position(16.853, -5.219, 1.725);
        atom_c1p.occ = 1.0f;
        atom_c1p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c1p);

        gemmi::Atom atom_c2p;
        atom_c2p.name = " C2'";
        atom_c2p.element = gemmi::Element("C");
        atom_c2p.pos = gemmi::Position(17.473, -4.888, 3.087);
        atom_c2p.occ = 1.0f;
        atom_c2p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c2p);

        gemmi::Atom atom_c3p;
        atom_c3p.name = " C3'";
        atom_c3p.element = gemmi::Element("C");
        atom_c3p.pos = gemmi::Position(16.678, -5.807, 4.033);
        atom_c3p.occ = 1.0f;
        atom_c3p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c3p);

        EXPECT_EQ(res_copy.atoms.size(), 3U);

        coot::ideal_rna::add_o2_prime_gemmi(&res_copy);

        EXPECT_EQ(res_copy.atoms.size(), 4U);

        bool has_o2p = false;
        for (const auto& atom : res_copy.atoms) {
            if (atom.name.find("O2'") != std::string::npos ||
                atom.name.find("O2*") != std::string::npos) {
                has_o2p = true;
                // Expected: o2p(8.870, -5.158, 3.018) + shift(10, 0, 0) = (18.870, -5.158, 3.018)
                EXPECT_NEAR(atom.pos.x, 18.87, 0.01);
                EXPECT_NEAR(atom.pos.y, -5.158, 0.01);
                EXPECT_NEAR(atom.pos.z, 3.018, 0.01);
            }
        }
        EXPECT_TRUE(has_o2p);
    }

    // Test case 3: Residue without all required atoms (should not add O2*)
    {
        gemmi::Residue res_copy;
        res_copy.name = "G";
        res_copy.seqid.num = 3;
        // res_copy.seqid.icode = ' ';

        // Only add C1' and C3' (missing C2')
        gemmi::Atom atom_c1p;
        atom_c1p.name = " C1'";
        atom_c1p.element = gemmi::Element("C");
        atom_c1p.pos = gemmi::Position(6.853, -5.219, 1.725);
        atom_c1p.occ = 1.0f;
        atom_c1p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c1p);

        gemmi::Atom atom_c3p;
        atom_c3p.name = " C3'";
        atom_c3p.element = gemmi::Element("C");
        atom_c3p.pos = gemmi::Position(6.678, -5.807, 4.033);
        atom_c3p.occ = 1.0f;
        atom_c3p.b_iso = 30.0f;
        res_copy.atoms.push_back(atom_c3p);

        EXPECT_EQ(res_copy.atoms.size(), 2U);

        coot::ideal_rna::add_o2_prime_gemmi(&res_copy);

        EXPECT_EQ(res_copy.atoms.size(), 2U);

        bool has_o2p = false;
        for (const auto& atom : res_copy.atoms) {
            if (atom.name.find("O2'") != std::string::npos ||
                atom.name.find("O2*") != std::string::npos) {
                has_o2p = true;
            }
        }
        EXPECT_FALSE(has_o2p);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}