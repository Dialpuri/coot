#pragma once

#include <gemmi/model.hpp>

namespace coot {

class molecule_t {
public:
    void change_chain_id_with_residue_range_helper_insert_or_add_gemmi(
        gemmi::Chain& chain, gemmi::Residue new_residue);
};

} // namespace coot

// Out-of-line definition (does not use `this` — operates only on parameters)
inline void coot::molecule_t::change_chain_id_with_residue_range_helper_insert_or_add_gemmi(
    gemmi::Chain& chain, gemmi::Residue new_residue)
{
    // Insert new_residue before the first residue whose seqnum > new_residue's seqnum.
    // If no such residue exists, append to the end.

    int resno_new_residue = new_residue.seqid.num.value;
    int target_index = -1;  // equivalent of RESIDUE_NUMBER_UNSET
    int best_seq_num_diff = 99999999;

    for (int i = 0; i < static_cast<int>(chain.residues.size()); i++) {
        int chain_residue_seq_num = chain.residues[i].seqid.num.value;
        int this_seq_num_diff = chain_residue_seq_num - resno_new_residue;
        if (this_seq_num_diff > 0) {
            if (this_seq_num_diff < best_seq_num_diff) {
                best_seq_num_diff = this_seq_num_diff;
                target_index = i;
            }
        }
    }

    if (target_index != -1) {
        // Insert before the found position
        chain.residues.insert(chain.residues.begin() + target_index, new_residue);
    } else {
        // Append to end
        chain.residues.push_back(new_residue);
    }
}