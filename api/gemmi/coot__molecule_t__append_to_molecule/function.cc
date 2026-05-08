#include "function.hh"
#include <string>
#include <vector>

namespace coot {

int append_to_molecule_gemmi(gemmi::Structure& dest, const gemmi::Structure& water_mol) {
    int istat = 0; // fail status initially.
    int n_atom = 0;  // 0 new atoms added initially.
    float default_new_atoms_b_factor = 20.0;

    // Get the first model from dest (gemmi models are 0-indexed)
    if (dest.models.empty()) {
        return istat;
    }
    
    gemmi::Model& dest_model = dest.models[0];
    
    // Run over the chains in water_mol
    for (const gemmi::Model& src_model : water_mol.models) {
        for (const gemmi::Chain& src_chain : src_model.chains) {
            std::string chain_id = src_chain.name;
            
            // Check if chain already exists in dest
            bool imatch = false;
            for (gemmi::Chain& dest_chain : dest_model.chains) {
                if (dest_chain.name == chain_id) {
                    imatch = true;
                    istat = 1;
                    break;
                }
            }
            
            if (!imatch) {
                // Create new chain
                gemmi::Chain new_chain;
                new_chain.name = chain_id;
                
                // Run over residues in the source chain
                for (const gemmi::Residue& src_res : src_chain.residues) {
                    if (src_res.atoms.empty()) {
                        continue;
                    }
                    
                    // Create new residue
                    gemmi::Residue new_res;
                    new_res.name = src_res.name;
                    new_res.seqid = src_res.seqid;
                    
                    // Add atoms
                    for (const gemmi::Atom& src_atom : src_res.atoms) {
                        gemmi::Atom new_atom;
                        new_atom.name = src_atom.name;
                        new_atom.element = src_atom.element;
                        new_atom.pos = src_atom.pos;
                        new_atom.occ = src_atom.occ;
                        new_atom.b_iso = default_new_atoms_b_factor;
                        new_res.atoms.push_back(new_atom);
                        n_atom++;
                    }
                    
                    new_chain.residues.push_back(new_res);
                }
                
                dest_model.chains.push_back(new_chain);
            }
        }
    }

    if (n_atom > 0) {
        // In gemmi, no FinishStructEdit() needed - just return
    }

    return istat;
}

}