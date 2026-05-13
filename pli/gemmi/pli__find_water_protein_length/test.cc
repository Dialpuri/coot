#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test suite for pli::find_water_protein_length
TEST(OracleTest, find_water_protein_length) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model* model = &st.models[0];

    // Case 1: Test with a protein residue (ILE at A/99) - should find nearby protein atoms
    {
        const gemmi::Chain* chain_a = nullptr;
        const gemmi::Residue* ligand_res = nullptr;
        for (const gemmi::Chain& chain : model->chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                for (const gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 99) {
                        ligand_res = &res;
                        break;
                    }
                }
                break;
            }
        }

        ASSERT_NE(ligand_res, nullptr) << "Residue ILE at A/99 not found";
        
        EXPECT_EQ(chain_a->name, "A");
        EXPECT_EQ(ligand_res->seqid.num.value, 99);
        EXPECT_EQ(ligand_res->name, "ILE");
        
        double dist = pli::find_water_protein_length_gemmi(ligand_res, model);
        EXPECT_NEAR(dist, 1.35638, 1e-4);
    }

    // Case 2: Test with water residue in chain B - should find protein atoms
    {
        const gemmi::Chain* chain_b = nullptr;
        const gemmi::Residue* water_res = nullptr;
        for (const gemmi::Chain& chain : model->chains) {
            if (chain.name == "B") {
                chain_b = &chain;
                for (const gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 1) {
                        water_res = &res;
                        break;
                    }
                }
                break;
            }
        }

        ASSERT_NE(water_res, nullptr) << "Residue HOH at B/1 not found";
        
        EXPECT_EQ(chain_b->name, "B");
        EXPECT_EQ(water_res->seqid.num.value, 1);
        EXPECT_EQ(water_res->name, "HOH");
        
        double dist = pli::find_water_protein_length_gemmi(water_res, model);
        EXPECT_NEAR(dist, 1.81198, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}