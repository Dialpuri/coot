#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MatchContainerForResiduesTest) {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
    ASSERT_GE(st.models.size(), 1) << "Failed to load PDB file or no models found";
    
    gemmi::Model& model = st.models[0];
    
    // Find chain A
    gemmi::Chain* chain = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr) << "Chain A not found";
    
    size_t n_res = chain->residues.size();
    ASSERT_GE(n_res, 11) << "Not enough residues in chain A";
    
    // Residue indices: GetResidue(9) -> index 9, GetResidue(10) -> index 10
    gemmi::Residue* res1 = &chain->residues[9];
    gemmi::Residue* res2 = &chain->residues[10];
    
    ASSERT_NE(res1, nullptr) << "Residue 10 not found";
    ASSERT_NE(res2, nullptr) << "Residue 11 not found";
    
    std::cout << "=== Test: Direct instantiation of match_container_for_residues_t_gemmi ===" << std::endl;
    std::cout << "INPUT res1: " << chain->name << "/" << res1->seqid.num.value << std::endl;
    std::cout << "INPUT res2: " << chain->name << "/" << res2->seqid.num.value << std::endl;
    
    coot::match_container_for_residues_t_gemmi container(res1, res2);
    
    // Test case 1: merge_flags.first = true (upstream merge)
    {
        std::cout << "=== Test case 1: merge_flags.first=true (upstream merge) ===" << std::endl;
        
        for (auto& at1 : res1->atoms) {
            std::string name1 = at1.name;
            
            for (auto& at2 : res2->atoms) {
                std::string name2 = at2.name;
                
                if (name1 == name2 && name1 != " O  ") {
                    container.add(&at1, &at2);
                }
            }
        }
        
        std::cout << "OUTPUT atom_pairs_size: " << container.atom_pairs.size() << std::endl;
        
        std::pair<bool, bool> merge_flags(true, false);
        
        std::cout << "INPUT merge_flags.first: " << merge_flags.first << std::endl;
        std::cout << "INPUT merge_flags.second: " << merge_flags.second << std::endl;
        
        size_t n_res_before = chain->residues.size();
        std::cout << "OUTPUT residues_before_meld: " << n_res_before << std::endl;
        container.meld_gemmi(model, merge_flags);
        size_t n_res_after = chain->residues.size();
        std::cout << "OUTPUT residues_after_meld: " << n_res_after << std::endl;
        std::cout << "OUTPUT residue_count_changed: " << (n_res_before != n_res_after ? "true" : "false") << std::endl;
        
        // Rebuild container for second test
        container = coot::match_container_for_residues_t_gemmi(res1, res2);
        for (auto& at1 : res1->atoms) {
            std::string name1 = at1.name;
            
            for (auto& at2 : res2->atoms) {
                std::string name2 = at2.name;
                
                if (name1 == name2 && name1 != " O  ") {
                    container.add(&at1, &at2);
                }
            }
        }
    }
    
    // Test case 2: merge_flags.first = false (downstream merge)
    {
        std::cout << "=== Test case 2: merge_flags.first=false (downstream merge) ===" << std::endl;
        
        std::pair<bool, bool> merge_flags(false, false);
        
        std::cout << "INPUT merge_flags.first: " << merge_flags.first << std::endl;
        std::cout << "INPUT merge_flags.second: " << merge_flags.second << std::endl;
        
        container.meld_gemmi(model, merge_flags);
        size_t n_res_after = chain->residues.size();
        std::cout << "OUTPUT residues_after_meld: " << n_res_after << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}