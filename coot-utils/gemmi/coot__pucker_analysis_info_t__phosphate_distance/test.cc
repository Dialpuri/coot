#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(PuckerAnalysisInfoTest, PhosphateDistance) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find nucleotide residues (not protein residues like MET or GLU)
    gemmi::Model& model = st.models[0];
    
    // Look for residue 1 (should be a nucleotide based on the test context)
    gemmi::Residue* res1 = nullptr;
    gemmi::Residue* res2 = nullptr;
    
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == 1) {
                res1 = &res;
            } else if (res.seqid.num.value == 2) {
                res2 = &res;
            }
        }
    }
    
    // Test case 1: Protein residue (MET) should throw
    if (res1 && res1->name == "MET") {
        char altconf = ' ';
        EXPECT_THROW({
            coot::pucker_analysis_info_t pai(*res1, altconf);
        }, std::exception);
    }
    
    // Test case 2: Protein residue (GLU) should throw
    if (res2 && res2->name == "GLU") {
        char altconf = ' ';
        EXPECT_THROW({
            coot::pucker_analysis_info_t pai(*res2, altconf);
        }, std::exception);
    }
    
    // Find a nucleotide residue to test phosphate_distance
    gemmi::Residue* nucleotide_res = nullptr;
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            // Check if it's a nucleotide by looking for P atom
            bool has_phosphate = false;
            for (const auto& atom : res.atoms) {
                if (atom.name == " P  ") {
                    has_phosphate = true;
                    break;
                }
            }
            if (has_phosphate) {
                nucleotide_res = &res;
                break;
            }
        }
        if (nucleotide_res) break;
    }
    
    if (nucleotide_res) {
        char altconf = ' ';
        coot::pucker_analysis_info_t pai(*nucleotide_res, altconf);
        
        // Test phosphate_distance - should work for nucleotide
        // The following residue would be the next one in the chain
        if (!nucleotide_res->atoms.empty()) {
            EXPECT_NO_THROW({
                float dist = pai.phosphate_distance(*nucleotide_res);
                // Just verify it doesn't throw - actual value depends on structure
                EXPECT_GE(dist, 0.0f);
            });
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}