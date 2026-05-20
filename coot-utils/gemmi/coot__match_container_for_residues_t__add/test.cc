#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, match_container_for_residues_add) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find residue A/5
    gemmi::Residue *res1 = nullptr;
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                if (res.seqid.num.value == 5 && chain.name == "A") {
                    res1 = &res;
                    break;
                }
    
    // Find residue A/6
    gemmi::Residue *res2 = nullptr;
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                if (res.seqid.num.value == 6 && chain.name == "A") {
                    res2 = &res;
                    break;
                }

    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    // Get CA atoms from each residue
    gemmi::Atom *at1 = nullptr;
    gemmi::Atom *at2 = nullptr;
    for (auto& atom : res1->atoms)
        if (atom.name == "CA") { at1 = &atom; break; }
    for (auto& atom : res2->atoms)
        if (atom.name == "CA") { at2 = &atom; break; }
    
    ASSERT_NE(at1, nullptr);
    ASSERT_NE(at2, nullptr);

    // Case 1: add a pair of CA atoms from A/5 and A/6
    {
        coot::match_container_for_residues_t m(res1, res2);
        EXPECT_EQ(m.atom_pairs.size(), 0);
        m.add_gemmi(at1, at2);
        EXPECT_EQ(m.atom_pairs.size(), 1);
        EXPECT_EQ(res1->seqid.num.value, 5);
        EXPECT_EQ(res2->seqid.num.value, 6);
    }

    // Case 2: add two atoms from the same residue (A/5)
    {
        gemmi::Atom *atN = nullptr;
        gemmi::Atom *atCA = nullptr;
        for (auto& atom : res1->atoms) {
            if (atom.name == "N") atN = &atom;
            if (atom.name == "CA") atCA = &atom;
        }
        ASSERT_NE(atN, nullptr);
        ASSERT_NE(atCA, nullptr);
        coot::match_container_for_residues_t m(res1, res1);
        EXPECT_EQ(m.atom_pairs.size(), 0);
        m.add_gemmi(atN, atCA);
        EXPECT_EQ(m.atom_pairs.size(), 1);
        EXPECT_EQ(res1->seqid.num.value, 5);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}