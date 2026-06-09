#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <map>
#include <iostream>

TEST(GemmiTest, make_h_bond_restraints_from_res_vec_auto) {
    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    
    // Get first model
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    
    // Gather all residue ids, names, and chain ids from the structure
    std::vector<gemmi::SeqId> seq_ids;
    std::vector<std::string> res_names;
    std::vector<std::string> chain_ids;
    
    // Original test only iterated over chain 0 (first chain)
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            // Original MMDB test iterated all residues (no atom filter)
            seq_ids.push_back(res.seqid);
            res_names.push_back(res.name);
            chain_ids.push_back(chain.name);
        }
        // Only use the first chain, matching original test
        break;
    }
    
    // Original test expected 267 residues
    EXPECT_EQ(seq_ids.size(), 267u);
    
    // Set up protein geometry
    coot::protein_geometry pg;
    pg.init_standard();
    
    // Call the ported function
    EXPECT_NO_THROW(
        coot::make_h_bond_restraints_from_res_vec_auto_gemmi(
            seq_ids, res_names, chain_ids, &st, pg
        )
    );
    
    // Print intermediate values for debugging
    std::cout << "DEBUG: seq_ids = " << seq_ids.size() << std::endl;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
