#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t_detail {

// Insert or add residue to chain, maintaining sequence order
// Returns the position where residue was inserted (or added at end)
// Uses gemmi::Chain::residues (vector<Residue>) which supports insert()
inline size_t change_chain_id_with_residue_range_helper_insert_or_add(
    gemmi::Chain& chain,
    const gemmi::Residue& new_residue) {
  
  const int resno_new_residue = new_residue.seqid.num.value;
  size_t best_seq_num_diff = 99999999;
  size_t target_pos = chain.residues.size();  // default: add at end
  
  for (size_t i = 0; i < chain.residues.size(); ++i) {
    int chain_residue_seq_num = chain.residues[i].seqid.num.value;
    int this_seq_num_diff = chain_residue_seq_num - resno_new_residue;
    if (this_seq_num_diff > 0) {
      if (static_cast<size_t>(this_seq_num_diff) < best_seq_num_diff) {
        best_seq_num_diff = static_cast<size_t>(this_seq_num_diff);
        target_pos = i;
      }
    }
  }
  
  if (target_pos < chain.residues.size()) {
    // Insert before target position
    chain.residues.insert(chain.residues.begin() + target_pos, new_residue);
  } else {
    // Add at end
    chain.residues.push_back(new_residue);
  }
  
  return target_pos;
}

} // namespace molecule_t_detail
} // namespace coot
