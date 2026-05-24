#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "coot-utils/merge-atom-selections.hh"

namespace coot {

// Gemmi-based adapter for match_container_for_residues_t
// This mirrors the MMDB version but uses gemmi types
struct match_container_for_residues_t_gemmi {
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    struct atom_pair { gemmi::Atom *at1, *at2; };
    std::vector<atom_pair> atom_pairs;
    
    match_container_for_residues_t_gemmi(gemmi::Residue *r1, gemmi::Residue *r2)
        : residue_1(r1), residue_2(r2) {}
    
    void add(gemmi::Atom *a1, gemmi::Atom *a2) {
        atom_pairs.push_back({a1, a2});
    }
    
    // Ported from MMDB version - uses gemmi APIs
    void meld_gemmi(gemmi::Model& model, std::pair<bool, bool> merge_flags);
    
    // Helper: get chain containing a residue
    static gemmi::Chain* get_chain_ptr(gemmi::Residue* res, gemmi::Model& model);
    
    // Helper: get contiguous residues from a starting residue
    static std::vector<gemmi::Residue*> residue_vector_from_residue(gemmi::Residue* res, gemmi::Chain& chain);
    
    // Helper: meld residues into target chain
    static void meld_residues_gemmi(std::vector<gemmi::Residue*>& res_vec,
                                     gemmi::Residue* ref_res,
                                     int res_no_delta,
                                     gemmi::Chain* to_chain_p,
                                     gemmi::Model& model);
};

} // namespace coot