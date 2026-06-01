#pragma once

#include <gemmi/model.hpp>
#include <cmath>
#include <string>
#include <cstdlib>

namespace coot {
namespace daca {

/*
 * Check if atom_cra.atom is a mainchain atom (N, CA, C, O) within
 * ±1 residue of ref_cra.residue in the same chain.
 *
 * gemmi equivalent of coot::daca::atom_is_neighbour_mainchain.
 */
inline bool atom_is_neighbour_mainchain_gemmi(
    const gemmi::CRA& atom_cra,
    const gemmi::CRA& ref_cra)
{
  bool status = false;

  if (atom_cra.chain && atom_cra.residue && atom_cra.atom &&
      ref_cra.chain && ref_cra.residue)
  {
    if (atom_cra.chain == ref_cra.chain)
    {
      std::ptrdiff_t idx_res_1 = ref_cra.residue - ref_cra.chain->residues.data();
      std::ptrdiff_t idx_res_2 = atom_cra.residue - atom_cra.chain->residues.data();
      std::ptrdiff_t idx_delta = std::abs(idx_res_2 - idx_res_1);

      if (idx_delta < 2)
      {
        std::string atom_name(atom_cra.atom->name);
        if (atom_name == "N")  { status = true; }
        if (atom_name == "CA") { status = true; }
        if (atom_name == "C")  { status = true; }
        if (atom_name == "O")  { status = true; }
      }
    }
  }

  return status;
}

} // namespace daca
} // namespace coot