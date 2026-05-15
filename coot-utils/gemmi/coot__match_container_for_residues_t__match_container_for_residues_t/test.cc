#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, match_container_for_residues_gemmi) {
    // Setup: load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "PDB load failed";
    gemmi::Model& model = st.models[0];

    // Find residue pointers by iterating chain A
    gemmi::Residue *res1 = nullptr;
    gemmi::Residue *res2 = nullptr;
    
    for (gemmi::Residue& res : model.find_chain("A")->residues) {
        if (res.seqid.num.value == 10 && res.name == "ILE") {
            res1 = &res;
        }
        if (res.seqid.num.value == 11 && res.name == "GLY") {
            res2 = &res;
        }
    }
    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    // === Case 1: default constructor (null residues) ===
    {
        coot::match_container_for_residues_t_gemmi m;
        EXPECT_EQ(m.residue_1, nullptr);
        EXPECT_EQ(m.residue_2, nullptr);
        EXPECT_EQ(m.atom_pairs.size(), 0);
        EXPECT_EQ(m.fragment_1_res_vec.size(), 0);
        EXPECT_EQ(m.fragment_2_res_vec.size(), 0);
    }

    // === Case 2: two-arg constructor with valid residues, then add ===
    {
        coot::match_container_for_residues_t_gemmi m(res1, res2);

        // residue_1: ILE 10
        EXPECT_NE(m.residue_1, nullptr);
        if (m.residue_1) {
            EXPECT_STREQ(m.residue_1->name.c_str(), "ILE");
            EXPECT_EQ(m.residue_1->seqid.num.value, 10);
        }

        // residue_2: GLY 11
        EXPECT_NE(m.residue_2, nullptr);
        if (m.residue_2) {
            EXPECT_STREQ(m.residue_2->name.c_str(), "GLY");
            EXPECT_EQ(m.residue_2->seqid.num.value, 11);
        }

        // atom_pairs before add: 0
        EXPECT_EQ(m.atom_pairs.size(), 0);
        // fragment_1 before add: 0
        EXPECT_EQ(m.fragment_1_res_vec.size(), 0);

        // Add matching CA atoms
        gemmi::Atom *ca1 = nullptr;
        gemmi::Atom *ca2 = nullptr;
        
        for (gemmi::Atom& atom : res1->atoms) {
            if (atom.name == "CA") {
                ca1 = &atom;
                break;
            }
        }
        for (gemmi::Atom& atom : res2->atoms) {
            if (atom.name == "CA") {
                ca2 = &atom;
                break;
            }
        }
        ASSERT_NE(ca1, nullptr);
        ASSERT_NE(ca2, nullptr);
        m.add(ca1, ca2);

        // atom_pairs after add: 1
        EXPECT_EQ(m.atom_pairs.size(), 1);
        // fragment_1 after add: 0
        EXPECT_EQ(m.fragment_1_res_vec.size(), 0);
        // fragment_2 after add: 0
        EXPECT_EQ(m.fragment_2_res_vec.size(), 0);

        if (!m.atom_pairs.empty()) {
            auto &p = m.atom_pairs[0];
            gemmi::Atom *atom1 = p.first;
            gemmi::Atom *atom2 = p.second;

            // atom1_spec: ILE 10  CA 
            EXPECT_STREQ(res1->name.c_str(), "ILE");
            EXPECT_EQ(res1->seqid.num.value, 10);
            EXPECT_STREQ(atom1->name.c_str(), "CA");

            // atom2_spec: GLY 11  CA 
            EXPECT_STREQ(res2->name.c_str(), "GLY");
            EXPECT_EQ(res2->seqid.num.value, 11);
            EXPECT_STREQ(atom2->name.c_str(), "CA");
        }
    }

    // === Case 3: ctor with one null residue (edge case) ===
    {
        coot::match_container_for_residues_t_gemmi m(res1, nullptr);

        // residue_1: ILE 10
        EXPECT_NE(m.residue_1, nullptr);
        if (m.residue_1) {
            EXPECT_STREQ(m.residue_1->name.c_str(), "ILE");
            EXPECT_EQ(m.residue_1->seqid.num.value, 10);
        }

        // residue_2: NULL
        EXPECT_EQ(m.residue_2, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}