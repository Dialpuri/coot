#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, fix_up_residue_and_atom_names) {
    // Test case 1: DNA residue with old-style names (O1P, O2P, C5*)
    {
        gemmi::Residue res;
        res.name = "DA";
        res.seqid.num.value = 1;

        gemmi::Atom at1;
        at1.name = " O1P";
        at1.pos = gemmi::Position(0.0, 0.0, 0.0);
        at1.element = gemmi::Element("O");
        res.atoms.push_back(at1);

        gemmi::Atom at2;
        at2.name = " O2P";
        at2.pos = gemmi::Position(1.0, 0.0, 0.0);
        at2.element = gemmi::Element("O");
        res.atoms.push_back(at2);

        gemmi::Atom at3;
        at3.name = " C5*";
        at3.pos = gemmi::Position(2.0, 0.0, 0.0);
        at3.element = gemmi::Element("C");
        res.atoms.push_back(at3);

        coot::ideal_rna::fix_up_residue_and_atom_names_gemmi(&res, true);

        bool o1p_fixed = false, o2p_fixed = false;
        for (const auto& atom : res.atoms) {
            if (atom.name == " OP1") o1p_fixed = true;
            if (atom.name == " OP2") o2p_fixed = true;
        }
        EXPECT_TRUE(o1p_fixed);
        EXPECT_TRUE(o2p_fixed);
    }

    // Test case 2: RNA residue with asterisk in name (O2*)
    {
        gemmi::Residue res2;
        res2.name = "A";
        res2.seqid.num.value = 1;

        gemmi::Atom at;
        at.name = " O2*";
        at.pos = gemmi::Position(0.0, 0.0, 0.0);
        at.element = gemmi::Element("O");
        res2.atoms.push_back(at);

        coot::ideal_rna::fix_up_residue_and_atom_names_gemmi(&res2, false);

        bool asterisk_fixed = false;
        for (const auto& atom : res2.atoms) {
            if (atom.name == " O2'") asterisk_fixed = true;
        }
        EXPECT_TRUE(asterisk_fixed);
    }

    // Test case 3: DNA DT residue with C5M
    {
        gemmi::Residue res3;
        res3.name = "DT";
        res3.seqid.num.value = 1;

        gemmi::Atom at;
        at.name = " C5M";
        at.pos = gemmi::Position(0.0, 0.0, 0.0);
        at.element = gemmi::Element("C");
        res3.atoms.push_back(at);

        coot::ideal_rna::fix_up_residue_and_atom_names_gemmi(&res3, true);

        bool c5m_fixed = false;
        for (const auto& atom : res3.atoms) {
            if (atom.name == " C7 ") c5m_fixed = true;
        }
        EXPECT_TRUE(c5m_fixed);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}