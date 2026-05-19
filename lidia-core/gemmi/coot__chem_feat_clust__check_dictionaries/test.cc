#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, check_dictionaries) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    
    coot::protein_geometry geom;
    
    // ── Case 1: Standard residues with dictionary entries (should pass) ──
    {
        std::vector<coot::chem_feat_clust::chem_feat_solvated_ligand> ligands;
        
        // Get residue 10 (A/10)
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& residue : chain.residues) {
                    if (residue.seqid.num.value == 10 && chain.name == "A") {
                        coot::chem_feat_clust::chem_feat_solvated_ligand cl;
                        cl.residue = const_cast<gemmi::Residue*>(&residue);
                        ligands.push_back(cl);
                        break;
                    }
                }
                if (!ligands.empty()) break;
            }
            if (!ligands.empty()) break;
        }
        
        // Get residue 11 (A/11)
        for (const auto& model : st.models) {
            for (const auto& chain : model.chains) {
                for (const auto& residue : chain.residues) {
                    if (residue.seqid.num.value == 11 && chain.name == "A") {
                        coot::chem_feat_clust::chem_feat_solvated_ligand cl;
                        cl.residue = const_cast<gemmi::Residue*>(&residue);
                        ligands.push_back(cl);
                        break;
                    }
                }
                if (ligands.size() == 2) break;
            }
            if (ligands.size() == 2) break;
        }
        
        bool result = coot::chem_feat_clust::check_dictionaries_gemmi(ligands, geom, 0);
        EXPECT_EQ(result, true);
    }
    
    // ── Case 2: Ligand with null residue pointer (guard clause path) ──
    {
        std::vector<coot::chem_feat_clust::chem_feat_solvated_ligand> ligands;
        
        coot::chem_feat_clust::chem_feat_solvated_ligand cl2;
        cl2.residue = nullptr;
        ligands.push_back(cl2);
        
        bool result = coot::chem_feat_clust::check_dictionaries_gemmi(ligands, geom, 0);
        EXPECT_EQ(result, true);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}