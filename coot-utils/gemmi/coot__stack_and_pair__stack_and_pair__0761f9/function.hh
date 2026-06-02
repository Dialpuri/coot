#pragma once

#include <map>
#include <utility>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__calculate_residue_normals__bc3097/gemmi/function.hh"

namespace coot {
namespace stack_and_pair {

// gemmi port of coot::stack_and_pair::stack_and_pair(mmdb::Manager*, int)
// The original constructor:
//   1. calls init() — resets angle_crit and normal_map
//   2. gets selected atoms from the manager by selection handle
//   3. calls calculate_residue_normals(sel_atoms, n_sel_atoms)
//
// In gemmi, there is no "selection handle" concept. We accept the whole
// model, collect all residues, and call the existing calculate_residue_normals_gemmi.
std::map<const gemmi::Residue*, clipper::Coord_orth>
stack_and_pair_gemmi(gemmi::Model& model) {
   // init() equivalent: start with empty map (no persistent state in this free function)

   // Collect all residues from the model (equivalent of GetSelIndex for all atoms)
   std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
   for (gemmi::Chain& chain : model.chains) {
      for (gemmi::Residue& residue : chain.residues) {
         residues_vec.push_back(std::make_pair(false, &residue));
      }
   }

   return calculate_residue_normals_gemmi(residues_vec);
}

} // namespace stack_and_pair
} // namespace coot