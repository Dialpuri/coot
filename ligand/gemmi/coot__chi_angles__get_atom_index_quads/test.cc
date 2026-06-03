#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chi_angles) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1u);

        gemmi::Model &model = st.models[0];

        // Get chains
        EXPECT_EQ(model.chains.size(), 2u);

        // Select chain A
        gemmi::Chain *chainA = nullptr;
        for (gemmi::Chain &ch : model.chains) {
            if (ch.name == "A") {
                chainA = &ch;
                break;
            }
        }
        ASSERT_NE(chainA, nullptr);

        // Get residues from chain A
        EXPECT_EQ(chainA->residues.size(), 267u);

        // Find a LYS residue (e.g., residue 24)
        gemmi::Residue *target_res = nullptr;
        for (gemmi::Residue &res : chainA->residues) {
            if (res.seqid.num.value == 24 && res.name == "LYS") {
                target_res = &res;
                break;
            }
        }

        if (!target_res) {
            // Fallback: find any LYS
            for (gemmi::Residue &res : chainA->residues) {
                if (res.name == "LYS") {
                    target_res = &res;
                    break;
                }
            }
        }

        ASSERT_NE(target_res, nullptr);

        std::string resname = target_res->name;
        int seqnum = target_res->seqid.num.value;
        EXPECT_EQ(resname, "LYS");
        EXPECT_EQ(seqnum, 24);

        // Get atoms from this residue
        EXPECT_EQ(target_res->atoms.size(), 9u);

        // Build atom name vector
        std::vector<std::string> residue_atom_names;
        for (gemmi::Atom &a : target_res->atoms) {
            residue_atom_names.push_back(a.name);
        }

        // Create atom_name_quads for LYS (hardcoded — no gemmi port for atom_name_quad_list)
        std::vector<coot::atom_name_quad> atom_name_quads;
        atom_name_quads.emplace_back(" N  ", " CA ", " CB ", " CG ");
        atom_name_quads.emplace_back(" CA ", " CB ", " CG ", " CD ");
        atom_name_quads.emplace_back(" CB ", " CG ", " CD ", " CE ");
        atom_name_quads.emplace_back(" CG ", " CD ", " CE ", " NZ ");

        EXPECT_EQ(atom_name_quads.size(), 4u);

        EXPECT_EQ(atom_name_quads[0].atom_name(0), " N  ");
        EXPECT_EQ(atom_name_quads[0].atom_name(1), " CA ");
        EXPECT_EQ(atom_name_quads[0].atom_name(2), " CB ");
        EXPECT_EQ(atom_name_quads[0].atom_name(3), " CG ");

        EXPECT_EQ(atom_name_quads[1].atom_name(0), " CA ");
        EXPECT_EQ(atom_name_quads[1].atom_name(1), " CB ");
        EXPECT_EQ(atom_name_quads[1].atom_name(2), " CG ");
        EXPECT_EQ(atom_name_quads[1].atom_name(3), " CD ");

        EXPECT_EQ(atom_name_quads[2].atom_name(0), " CB ");
        EXPECT_EQ(atom_name_quads[2].atom_name(1), " CG ");
        EXPECT_EQ(atom_name_quads[2].atom_name(2), " CD ");
        EXPECT_EQ(atom_name_quads[2].atom_name(3), " CE ");

        EXPECT_EQ(atom_name_quads[3].atom_name(0), " CG ");
        EXPECT_EQ(atom_name_quads[3].atom_name(1), " CD ");
        EXPECT_EQ(atom_name_quads[3].atom_name(2), " CE ");
        EXPECT_EQ(atom_name_quads[3].atom_name(3), " NZ ");

        // Call get_atom_index_quads_gemmi
        std::vector<coot::atom_index_quad> index_quads =
            coot::chi_angles::get_atom_index_quads_gemmi(atom_name_quads, residue_atom_names);

        EXPECT_EQ(index_quads.size(), 4u);

        EXPECT_EQ(index_quads[0].index1, 0);
        EXPECT_EQ(index_quads[0].index2, 1);
        EXPECT_EQ(index_quads[0].index3, 4);
        EXPECT_EQ(index_quads[0].index4, 5);

        EXPECT_EQ(index_quads[1].index1, 1);
        EXPECT_EQ(index_quads[1].index2, 4);
        EXPECT_EQ(index_quads[1].index3, 5);
        EXPECT_EQ(index_quads[1].index4, 6);

        EXPECT_EQ(index_quads[2].index1, 4);
        EXPECT_EQ(index_quads[2].index2, 5);
        EXPECT_EQ(index_quads[2].index3, 6);
        EXPECT_EQ(index_quads[2].index4, 7);

        EXPECT_EQ(index_quads[3].index1, 5);
        EXPECT_EQ(index_quads[3].index2, 6);
        EXPECT_EQ(index_quads[3].index3, 7);
        EXPECT_EQ(index_quads[3].index4, 8);

        // Edge case: empty atom_name_quads
        {
            std::vector<coot::atom_name_quad> empty_quads;
            std::vector<coot::atom_index_quad> empty_result =
                coot::chi_angles::get_atom_index_quads_gemmi(empty_quads, residue_atom_names);
            EXPECT_EQ(empty_result.size(), 0u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}