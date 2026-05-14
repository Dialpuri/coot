#include "function.hh"
#include <gemmi/pdb.hpp>

namespace coot {

void Bond_lines_container::atom_selection_missing_loops_gemmi(const atom_selection_container_t &asc,
                                                               int udd_atom_index_handle,
                                                               int udd_fixed_during_refinement_handle) {
    if (!asc.model) return;
    
    for (gemmi::Chain& chain : asc.model->chains) {
        int nres = chain.residues.size();
        if (nres < 2) continue;
        
        for (size_t ires = 1; ires < static_cast<size_t>(nres); ires++) {
            gemmi::Residue& residue_prev = chain.residues[ires - 1];
            gemmi::Residue& residue_this = chain.residues[ires];
            
            int n_atoms_prev = residue_prev.atoms.size();
            int n_atoms_this = residue_this.atoms.size();
            
            if (n_atoms_prev == 0 || n_atoms_this == 0) continue;
            
            int res_no_1 = residue_prev.seqid.num.value;
            int res_no_2 = residue_this.seqid.num.value;
            int res_no_delta = res_no_2 - res_no_1;
            
            if (res_no_delta > 1) {
                // Gap found - would call do_Ca_loop here with gemmi parameters
            }
        }
    }
}

} // namespace coot