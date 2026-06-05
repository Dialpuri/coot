#include "function.hh"
#include <gemmi/calculate.hpp>

namespace coot {
namespace molecule_t {

int full_atom_spec_to_atom_index_gemmi(
    const coot::atom_spec_t &atom_spec,
    gemmi::Structure &st) {

  // Normalize insertion code: MMDB uses "" for "none", gemmi uses ' '
  auto norm = [](const std::string &ic) { return ic.empty() ? std::string(" ") : ic; };
  std::string norm_ins_code = norm(atom_spec.ins_code);

  if (st.models.empty()) return -1;
  gemmi::Model &model = st.models[0];

  int idx = 0;
  for (gemmi::Chain &chain : model.chains) {
    if (chain.name != atom_spec.chain_id) {
      idx += static_cast<int>(gemmi::count_atom_sites(chain));
      continue;
    }
    for (gemmi::Residue &residue : chain.residues) {
      if (residue.seqid.num.value != atom_spec.res_no) {
        idx += static_cast<int>(gemmi::count_atom_sites(residue));
        continue;
      }
      if (std::string(1, residue.seqid.icode) != norm_ins_code) {
        idx += static_cast<int>(gemmi::count_atom_sites(residue));
        continue;
      }
      for (gemmi::Atom &atom : residue.atoms) {
        // altloc comparison
        char alt = atom.altloc;
        if (!atom_spec.alt_conf.empty()) {
          if (alt != atom_spec.alt_conf[0]) {
            idx++;
            continue;
          }
        } else {
          // When alt_conf is empty, match only atoms with no altloc
          // (gemmi uses '\0' for no altloc, but PDB files typically use ' ')
          // Follow coot convention: empty alt_conf matches altloc == ' ' or '\0'
          if (alt != ' ' && alt != '\0') {
            idx++;
            continue;
          }
        }
        if (atom.name == atom_spec.atom_name) {
          return idx;
        }
        idx++;
      }
    }
  }
  return -1;
}

} // namespace molecule_t
} // namespace coot