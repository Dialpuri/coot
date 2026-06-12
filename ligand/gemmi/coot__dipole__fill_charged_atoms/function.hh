#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__dipole__charged_atoms__df4eac/gemmi/function.hh"

namespace coot {
namespace dipole {

void fill_charged_atoms_gemmi(
    gemmi::Residue *residue_p,
    const coot::dictionary_residue_restraints_t &rest) {

  std::vector<std::pair<gemmi::Atom *, float>> v =
      charged_atoms_gemmi(residue_p, rest);

  for (unsigned int i = 0; i < v.size(); i++) {
    v[i].first->charge = static_cast<int>(v[i].second * 100.0f);
  }
}

} // namespace dipole
} // namespace coot
