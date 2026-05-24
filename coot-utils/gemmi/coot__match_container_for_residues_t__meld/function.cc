#include "function.hh"

namespace coot {

void
match_container_for_residues_t_gemmi::meld_gemmi(gemmi::Model& model, std::pair<bool, bool> merge_flags) {
    
    if (true) {
        // merge 2 into 1
        if (merge_flags.first) {
            // Get chain pointers - need to find the chain containing each residue
            gemmi::Chain* to_chain_p = get_chain_ptr(residue_1, model);
            if (!to_chain_p) return;
            
            int res_no_delta = residue_1->seqid.num.value - residue_2->seqid.num.value;
            
            // Get contiguous residue vector from residue_2
            gemmi::Chain* chain2 = get_chain_ptr(residue_2, model);
            if (!chain2) return;
            std::vector<gemmi::Residue*> res_vec = residue_vector_from_residue(residue_2, *chain2);
            
            // Meld the residues
            meld_residues_gemmi(res_vec, residue_2, res_no_delta, to_chain_p, model);
            
        } else {
            // downstream merge
            gemmi::Chain* to_chain_p = get_chain_ptr(residue_2, model);
            if (!to_chain_p) return;
            
            int res_no_delta = residue_2->seqid.num.value - residue_1->seqid.num.value;
            
            // Get contiguous residue vector from residue_1
            gemmi::Chain* chain1 = get_chain_ptr(residue_1, model);
            if (!chain1) return;
            std::vector<gemmi::Residue*> res_vec = residue_vector_from_residue(residue_1, *chain1);
            
            // Adjust residue numbers
            for (auto* res : res_vec) {
                res->seqid.num.value += res_no_delta;
            }
            
            // Meld the residues
            meld_residues_gemmi(res_vec, residue_1, 0, to_chain_p, model);
            
            // Ensure first residue number >= 1
            int n_residues = to_chain_p->residues.size();
            if (n_residues > 0) {
                int res_no_first = to_chain_p->residues[0].seqid.num.value;
                if (res_no_first < 1) {
                    int res_no_delta = 1 - res_no_first;
                    for (auto& res : to_chain_p->residues) {
                        res.seqid.num.value += res_no_delta;
                    }
                }
            }
        }
    }
}

// Helper: find chain containing a residue
// In gemmi, Residue has no parent pointer, so we search through all chains
gemmi::Chain*
match_container_for_residues_t_gemmi::get_chain_ptr(gemmi::Residue* res, gemmi::Model& model) {
    for (auto& chain : model.chains) {
        for (auto& r : chain.residues) {
            if (&r == res) {
                return &chain;
            }
        }
    }
    return nullptr;
}

// Helper: get contiguous residues from a starting residue
// This mimics MMDB's behavior - finds residues in sequence order
std::vector<gemmi::Residue*>
match_container_for_residues_t_gemmi::residue_vector_from_residue(gemmi::Residue* res, gemmi::Chain& chain) {
    std::vector<gemmi::Residue*> result;
    
    // Find the index of res in chain.residues
    size_t start_idx = 0;
    bool found = false;
    for (size_t i = 0; i < chain.residues.size(); i++) {
        if (&chain.residues[i] == res) {
            start_idx = i;
            found = true;
            break;
        }
    }
    
    if (!found) return result;
    
    // Add all residues from start_idx onwards
    for (size_t i = start_idx; i < chain.residues.size(); i++) {
        result.push_back(&chain.residues[i]);
    }
    
    return result;
}

// Helper: meld residues into target chain
void
match_container_for_residues_t_gemmi::meld_residues_gemmi(std::vector<gemmi::Residue*>& res_vec,
                                                     gemmi::Residue* ref_res,
                                                     int res_no_delta,
                                                     gemmi::Chain* to_chain_p,
                                                     gemmi::Model& model) {
    for (auto* res : res_vec) {
        // Calculate new residue number
        int new_seq_num = res->seqid.num.value + res_no_delta;
        
        // Find if residue with this number already exists in target chain
        bool found = false;
        size_t insert_idx = to_chain_p->residues.size();
        for (size_t i = 0; i < to_chain_p->residues.size(); i++) {
            if (to_chain_p->residues[i].seqid.num.value == new_seq_num) {
                // Merge atoms into existing residue
                for (auto& atom : res->atoms) {
                    // Deep copy atom
                    gemmi::Atom new_atom = atom;
                    new_atom.serial = static_cast<int>(to_chain_p->residues.back().atoms.size() + 1);
                    to_chain_p->residues[i].atoms.push_back(new_atom);
                }
                found = true;
                break;
            } else if (to_chain_p->residues[i].seqid.num.value > new_seq_num) {
                insert_idx = i;
                break;
            }
        }
        
        if (!found) {
            // Insert at correct position or append
            gemmi::Residue new_res = *res;
            new_res.seqid.num.value = new_seq_num;
            
            // Update atom serials
            int serial_offset = 0;
            if (!to_chain_p->residues.empty()) {
                for (const auto& r : to_chain_p->residues) {
                    for (const auto& a : r.atoms) {
                        serial_offset = std::max(serial_offset, a.serial);
                    }
                }
            }
            for (auto& atom : new_res.atoms) {
                atom.serial = ++serial_offset;
            }
            
            // Insert at correct position
            if (insert_idx < to_chain_p->residues.size()) {
                to_chain_p->residues.insert(to_chain_p->residues.begin() + insert_idx, new_res);
            } else {
                to_chain_p->residues.push_back(new_res);
            }
        }
    }
}

} // namespace coot