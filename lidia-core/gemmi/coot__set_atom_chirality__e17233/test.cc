#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, set_atom_chirality) {
    // Create an RDKit carbon atom (atomic number 6) — a typical chiral centre
    RDKit::Atom *rdkit_at = new RDKit::Atom();
    rdkit_at->setAtomicNum(6);

    // --- Case 1: pdbx_stereo_config = "R" ---
    {
        coot::dict_atom dict_atom;
        dict_atom.atom_id = "CA";
        dict_atom.atom_id_4c = " CA ";
        dict_atom.pdbx_stereo_config = std::make_pair(true, "R");

        RDKit::Atom::ChiralType tag_before = rdkit_at->getChiralTag();
        std::string cip_before;
        try { rdkit_at->getProp("_CIPCode", cip_before); } catch (...) { cip_before = "(none)"; }

        coot::set_atom_chirality_gemmi(rdkit_at, dict_atom);

        RDKit::Atom::ChiralType tag_after = rdkit_at->getChiralTag();
        std::string cip_after;
        try { rdkit_at->getProp("_CIPCode", cip_after); } catch (...) { cip_after = "(none)"; }

        EXPECT_EQ(static_cast<int>(tag_before), 0);
        EXPECT_EQ(static_cast<int>(tag_after), 1);
        EXPECT_EQ(cip_after, "R");
    }

    // Reset atom
    rdkit_at->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED);

    // --- Case 2: pdbx_stereo_config = "S" ---
    {
        coot::dict_atom dict_atom;
        dict_atom.atom_id = "CA";
        dict_atom.atom_id_4c = " CA ";
        dict_atom.pdbx_stereo_config = std::make_pair(true, "S");

        RDKit::Atom::ChiralType tag_before = rdkit_at->getChiralTag();

        coot::set_atom_chirality_gemmi(rdkit_at, dict_atom);

        RDKit::Atom::ChiralType tag_after = rdkit_at->getChiralTag();
        std::string cip_after;
        try { rdkit_at->getProp("_CIPCode", cip_after); } catch (...) { cip_after = "(none)"; }

        EXPECT_EQ(static_cast<int>(tag_before), 0);
        EXPECT_EQ(static_cast<int>(tag_after), 2);
        EXPECT_EQ(cip_after, "S");
    }

    // Reset atom
    rdkit_at->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED);

    // --- Case 3: pdbx_stereo_config flag = false (no stereo info) ---
    {
        coot::dict_atom dict_atom;
        dict_atom.atom_id = "CA";
        dict_atom.atom_id_4c = " CA ";
        dict_atom.pdbx_stereo_config = std::make_pair(false, "");

        RDKit::Atom::ChiralType tag_before = rdkit_at->getChiralTag();

        coot::set_atom_chirality_gemmi(rdkit_at, dict_atom);

        RDKit::Atom::ChiralType tag_after = rdkit_at->getChiralTag();

        EXPECT_EQ(static_cast<int>(tag_before), 0);
        EXPECT_EQ(static_cast<int>(tag_after), 0);
    }

    delete rdkit_at;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
