#pragma once

#include <gemmi/model.hpp>

namespace coot {

// Port of MMDB's delete_the_matched_residues_matched_residue to gemmi.
// Deletes either residue_1 or residue_2 from the structure based on the flag.
// short_fragment_is_in_first_selection = false means delete residue_1
// short_fragment_is_in_first_selection = true means delete residue_2
void delete_the_matched_residues_matched_residue_gemmi(gemmi::Structure& st,
                                                       gemmi::Residue* residue_1,
                                                       gemmi::Residue* residue_2,
                                                       bool short_fragment_is_in_first_selection) {
    // Determine which residue to delete
    gemmi::Residue* residue_to_delete = short_fragment_is_in_first_selection ? residue_2 : residue_1;
    
    // Find and delete the residue from its chain
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (size_t i = 0; i < chain.residues.size(); ++i) {
                if (&chain.residues[i] == residue_to_delete) {
                    chain.residues.erase(chain.residues.begin() + i);
                    break;  // Only delete one residue and exit
                }
            }
        }
    }
}

}