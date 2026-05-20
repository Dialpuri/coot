#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/contact.hpp>
#include "function.hh"

TEST(OracleTest, delete_downstream_selections) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get first model
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    
    // Build atom vectors for two selections:
    // - Selection 1: residues 10-12 (upstream of overlap)
    // - Selection 2: residues 11-17 (overlaps and extends downstream)
    std::vector<gemmi::CRA> atoms_sel1;
    std::vector<gemmi::CRA> atoms_sel2;
    
    for (auto& chain : model.chains) {
        if (chain.name != "A") continue;
        
        for (auto& res : chain.residues) {
            int seqnum = res.seqid.num.value;
            gemmi::CRA cra{&chain, &res, nullptr};
            
            if (seqnum >= 10 && seqnum <= 12) {
                for (auto& atom : res.atoms) {
                    cra.atom = &atom;
                    atoms_sel1.push_back(cra);
                }
            }
            if (seqnum >= 11 && seqnum <= 17) {
                for (auto& atom : res.atoms) {
                    cra.atom = &atom;
                    atoms_sel2.push_back(cra);
                }
            }
        }
    }
    
    // Build atom_pairs with matching atom names - use the exact same logic as MMDB test
    coot::match_container_for_residues_t mcr;
    double max_dist = 0.8;
    
    // Use gemmi's neighbor search
    gemmi::NeighborSearch ns(model, st.cell, 5.0);
    ns.populate(false);
    
    // Use gemmi's contact search to find contacts
    gemmi::ContactSearch cs(max_dist);
    cs.ignore = gemmi::ContactSearch::Ignore::Nothing;
    cs.setup_atomic_radii(1.0, 0.0);  // VdW radii
    
    // Find all contacts and filter
    std::vector<gemmi::ContactSearch::Result> all_contacts = cs.find_contacts(ns);
    
    for (const auto& c : all_contacts) {
        // Check if partner1 is in sel1 and partner2 is in sel2
        bool p1_in_sel1 = false;
        bool p2_in_sel2 = false;
        
        for (auto& cra : atoms_sel1) {
            if (cra.chain == c.partner1.chain && 
                cra.residue == c.partner1.residue && 
                cra.atom == c.partner1.atom) {
                p1_in_sel1 = true;
                break;
            }
        }
        
        for (auto& cra : atoms_sel2) {
            if (cra.chain == c.partner2.chain && 
                cra.residue == c.partner2.residue && 
                cra.atom == c.partner2.atom) {
                p2_in_sel2 = true;
                break;
            }
        }
        
        if (p1_in_sel1 && p2_in_sel2) {
            if (c.partner1.atom && c.partner2.atom) {
                std::string atom_name_1 = c.partner1.atom->name;
                std::string atom_name_2 = c.partner2.atom->name;
                
                if (atom_name_1 == atom_name_2) {
                    if (atom_name_1 != " O  ") {  // Skip oxygens
                        mcr.add(c.partner1, c.partner2);
                    }
                }
            }
        }
    }
    
    // Get initial residue count
    int res_count_before = 0;
    for (const auto& chain : model.chains) {
        res_count_before += chain.residues.size();
    }
    
    // Call delete_downstream with from_first=false (delete downstream from second selection)
    bool from_first = false;
    mcr.delete_downstream_gemmi(st, from_first, atoms_sel1, atoms_sel2);
    
    // Get residue count after
    int res_count_after = 0;
    for (const auto& chain : model.chains) {
        res_count_after += chain.residues.size();
    }
    
    // Assertions based on oracle output
    EXPECT_EQ(res_count_before, 268);
    EXPECT_EQ(static_cast<int>(mcr.atom_pairs.size()), 0);
    EXPECT_EQ(res_count_after, 268);
    EXPECT_EQ(res_count_before - res_count_after, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}