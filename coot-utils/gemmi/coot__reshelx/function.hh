#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/calculate.hpp>
#include <string>
#include "coot/coot-utils/coot-shelx.hh"

namespace coot {

// Port of coot::reshelx from MMDB to gemmi
// Takes a gemmi::Structure and returns a new Structure with reshelx transformations
inline gemmi::Structure reshelx_gemmi(const gemmi::Structure& st) {
    gemmi::Structure shelx_st;
    
    // Create a model and chain in the new structure
    gemmi::Model model;
    gemmi::Chain chain("A");
    
    bool made_afix_transfer_message = false;
    
    // Run over chains of the existing structure
    if (!st.models.empty()) {
        const gemmi::Model& model_p = st.models[0];
        int nchains_local = static_cast<int>(model_p.chains.size());
        
        for (int ichain = 0; ichain < nchains_local; ichain++) {
            const gemmi::Chain& chain_p = model_p.chains[ichain];
            int nres = static_cast<int>(chain_p.residues.size());
            
            for (int ires = 0; ires < nres; ires++) {
                const gemmi::Residue& residue_p = chain_p.residues[ires];
                
                // Deep copy the residue
                gemmi::Residue copy_residue = residue_p;
                
                // Copy atoms
                copy_residue.atoms = residue_p.atoms;
                
                chain.residues.push_back(copy_residue);
            }
        }
    }
    
    model.chains.push_back(chain);
    shelx_st.models.push_back(model);
    
    // Cell info - set all parameters to 0 to indicate no cell
    shelx_st.cell.a = 0;
    shelx_st.cell.b = 0;
    shelx_st.cell.c = 0;
    shelx_st.cell.alpha = 0;
    shelx_st.cell.beta = 0;
    shelx_st.cell.gamma = 0;
    
    // Copy spacegroup
    if (!st.spacegroup_hm.empty()) {
        shelx_st.spacegroup_hm = st.spacegroup_hm;
    }
    
    return shelx_st;
}

} // namespace coot