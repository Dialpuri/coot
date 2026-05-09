#include "function.hh"

namespace coot {

// Helper function to find atom in gemmi model by spec
inline gemmi::Atom* find_atom(gemmi::Model& model, const residue_spec_t& spec, const std::string& atom_name) {
    for (auto& chain : model.chains) {
        if (chain.name == spec.chain_name) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == spec.seq_num && res.seqid.icode == spec.ins_code) {
                    for (auto& atom : res.atoms) {
                        if (atom.name == atom_name) {
                            return &atom;
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

// Implementation of replace_fragment_gemmi
int replace_fragment_gemmi(gemmi::Model& model, atom_selection_container_t asc) {
    if (!asc.model) return 0;
    
    // Make backup of current state
    make_backup_gemmi("replace_fragment");
    
    for (const auto& atom_ctx : asc.atoms_with_context) {
        gemmi::Atom* at = std::get<0>(atom_ctx);
        std::string chain_name = std::get<1>(atom_ctx);
        int seq_num = std::get<2>(atom_ctx);
        char ins_code = std::get<3>(atom_ctx);
        
        // Try to find existing atom by spec
        int idx = -1;
        
        // Check if we can find by index
        if (asc.UDDOldAtomIndexHandle >= 0) {
            // In gemmi we don't have UDD handles, so skip this
        }
        
        if (idx == -1) {
            // Fall back to full spec search
            residue_spec_t spec(chain_name, seq_num, ins_code);
            gemmi::Atom* existing_atom = find_atom(model, spec, at->name);
            
            if (existing_atom) {
                // Update existing atom coordinates
                existing_atom->pos = at->pos;
                existing_atom->occ = at->occ;
                existing_atom->b_iso = at->b_iso;
                existing_atom->element = at->element;
                existing_atom->altloc = at->altloc;
            } else {
                // Add new atom - find or create chain and residue
                gemmi::Chain* chain_p = get_chain(model, chain_name);
                gemmi::Residue* residue_p = get_residue(model, residue_spec_t(chain_name, seq_num, ins_code));
                
                if (!chain_p) {
                    // Create new chain
                    chain_p = new gemmi::Chain();
                    chain_p->name = chain_name;
                    
                    // Create new residue
                    residue_p = new gemmi::Residue();
                    residue_p->name = at->name;  // Use the atom name as a placeholder for residue name
                    
                    // Create SeqId with proper value
                    residue_p->seqid.num = seq_num;
                    residue_p->seqid.icode = ins_code;
                    
                    chain_p->residues.push_back(*residue_p);
                    model.chains.push_back(*chain_p);
                } else {
                    if (!residue_p) {
                        // Create new residue
                        residue_p = new gemmi::Residue();
                        residue_p->name = at->name;  // Use the atom name as a placeholder for residue name
                        
                        // Create SeqId with proper value
                        residue_p->seqid.num = seq_num;
                        residue_p->seqid.icode = ins_code;
                        
                        int sn = find_serial_number_for_insert_gemmi(seq_num, std::string(1, ins_code), chain_name);
                        
                        if (sn != -1) {
                            // Insert at position - simplified
                            if (sn <= (int)chain_p->residues.size()) {
                                chain_p->residues.insert(chain_p->residues.begin() + sn, *residue_p);
                            } else {
                                chain_p->residues.push_back(*residue_p);
                            }
                        } else {
                            chain_p->residues.push_back(*residue_p);
                        }
                    }
                }
                
                if (residue_p) {
                    // Add atom to residue
                    gemmi::Atom new_atom;
                    new_atom.name = at->name;
                    new_atom.pos = at->pos;
                    new_atom.occ = at->occ;
                    new_atom.b_iso = at->b_iso;
                    new_atom.element = at->element;
                    new_atom.altloc = at->altloc;
                    
                    residue_p->atoms.push_back(new_atom);
                }
            }
        } else {
            // idx != -1 means we found an existing atom to update
            // In gemmi, we'd need to map index to atom somehow
            // For now, skip - this would require more complex indexing
        }
    }
    
    // Clean up - in real coot this would delete selection and reindex
    // For gemmi, we just return success
    
    return 1;
}

} // namespace coot