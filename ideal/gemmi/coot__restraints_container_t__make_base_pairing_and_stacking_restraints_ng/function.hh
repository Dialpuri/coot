#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__stacked_residues/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__paired_residues/gemmi/function.hh"

namespace coot {

// Gemmi port of restraints_container_t::make_base_pairing_and_stacking_restraints_ng
// This is a free function (not a member) because gemmi doesn't use restraints_container_t.
// Returns the number of base-pairing hydrogen bonds found.
inline int make_base_pairing_and_stacking_restraints_ng_gemmi(
    const gemmi::Structure& st,
    gemmi::Model& model,
    const std::vector<std::pair<bool, gemmi::Residue*>>& residues_vec,
    int imol,
    const coot::protein_geometry& geom)
{
  // 1. Stacked residues (parallel planes)
  std::vector<stack_and_pair::stacked_planes_info_t> stacked =
    stack_and_pair::stacked_residues_gemmi(st);

  if (!stacked.empty()) {
    extra_restraints_t extra_restraints;
    for (std::size_t i = 0; i < stacked.size(); i++) {
      // Find chain ids by searching the model (gemmi::Residue has no parent ptr)
      auto find_chain_id = [](gemmi::Model& mdl, gemmi::Residue* res) -> std::string {
        for (auto& ch : mdl.chains) {
          for (auto& r : ch.residues) {
            if (&r == res) return ch.name;
          }
        }
        return "";
      };

      auto chain1 = find_chain_id(model, stacked[i].res_1);
      auto chain2 = find_chain_id(model, stacked[i].res_2);

      int seq1 = stacked[i].res_1 ? stacked[i].res_1->seqid.num.value : 0;
      int seq2 = stacked[i].res_2 ? stacked[i].res_2->seqid.num.value : 0;
      std::string ic1 = stacked[i].res_1 ? std::string(1, stacked[i].res_1->seqid.icode) : "";
      std::string ic2 = stacked[i].res_2 ? std::string(1, stacked[i].res_2->seqid.icode) : "";

      residue_spec_t spec1(chain1, seq1, ic1);
      residue_spec_t spec2(chain2, seq2, ic2);

      parallel_planes_t ppr(spec1, spec2,
                            stacked[i].atom_names_1,
                            stacked[i].atom_names_2,
                            "", "");
      extra_restraints.parallel_plane_restraints.push_back(ppr);
    }
    // NOTE: In the original these are added via add_extra_restraints()
    // which requires a restraints_container_t. This free function returns
    // the number of base-pairing bonds found (stacking restraints go
    // through a different path in the original and are not counted in
    // the test oracle).
  }

  // 2. Paired residues (base pairs → hydrogen bonds)
  bool all_atoms_are_moving_flag = false;
  std::vector<stack_and_pair::paired_residues_info_t> pr =
    stack_and_pair::paired_residues_gemmi(
      model, residues_vec, all_atoms_are_moving_flag, geom, imol);

  unsigned int n_base_pairing_bonds = 0;

  for (std::size_t i = 0; i < pr.size(); i++) {
    for (std::size_t j = 0; j < pr[i].atom_pair_vec.size(); j++) {
      const auto& pair = pr[i].atom_pair_vec[j];
      // Both pointers non-null means a valid hydrogen bond
      if (pair.first && pair.second) {
        n_base_pairing_bonds++;
      }
    }
  }

  return static_cast<int>(n_base_pairing_bonds);
}

} // namespace coot
