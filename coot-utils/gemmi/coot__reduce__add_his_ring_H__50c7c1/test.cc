#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <string>

TEST(OracleTest, add_his_ring_H) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Valid HIS residue with second_neighb_vec size == 2
    {
        gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);

        gemmi::Residue* residue_p = chain->find_residue({gemmi::SeqId{60, ' '}, "", "HIS"});
        EXPECT_NE(residue_p, nullptr);

        EXPECT_EQ(residue_p->name, "HIS");
        EXPECT_EQ(residue_p->seqid.num.value, 60);

        int atoms_before = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_before, 10);

        std::string H_at_name = " HE2";
        std::string first_neigh = " CE1";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back("NE2");
        second_neighb_vec.push_back(" CD2");
        double bl = 0.86;

        coot::reduce::add_his_ring_H_gemmi(H_at_name, first_neigh, second_neighb_vec, bl, *residue_p);

        int atoms_after = static_cast<int>(residue_p->atoms.size());
        EXPECT_EQ(atoms_after, 11);

        // Check if HE2 atom was created
        bool atom_found = false;
        for (int i = 0; i < atoms_after; i++) {
            const auto& a = residue_p->atoms[i];
            std::string an = a.name;
            if (an == " HE2" || an == "HE2 ") {
                atom_found = true;
                EXPECT_EQ(an, " HE2");
                break;
            }
        }
        EXPECT_TRUE(atom_found);
    }

    // Case 2: second_neighb_vec.size() != 2 → guard clause, no-op
    {
        gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);

        gemmi::Residue* residue_p2 = chain->find_residue({gemmi::SeqId{71, ' '}, "", "HIS"});
        EXPECT_NE(residue_p2, nullptr);

        EXPECT_EQ(residue_p2->name, "HIS");
        EXPECT_EQ(residue_p2->seqid.num.value, 71);

        std::string H_at_name2 = " HE2";
        std::string first_neigh2 = " CE1";
        std::vector<std::string> second_neighb_vec2;
        second_neighb_vec2.push_back("NE2"); // size=1, guard clause skips

        int atoms_before2 = static_cast<int>(residue_p2->atoms.size());
        EXPECT_EQ(atoms_before2, 10);

        coot::reduce::add_his_ring_H_gemmi(H_at_name2, first_neigh2, second_neighb_vec2, 0.86, *residue_p2);

        int atoms_after2 = static_cast<int>(residue_p2->atoms.size());
        EXPECT_EQ(atoms_after2, 10);
        EXPECT_TRUE(atoms_before2 == atoms_after2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}