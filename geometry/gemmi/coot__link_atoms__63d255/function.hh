#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

inline std::pair<atom_spec_t, atom_spec_t>
link_atoms_gemmi(const gemmi::Connection& link, int model_number = 1) {

  const auto& a1_addr = link.partner1;
  const auto& a2_addr = link.partner2;

  // Convert gemmi icode (' ' = no insertion code) to coot convention ("")
  auto icode_to_str = [](char icode) -> std::string {
    return (icode == ' ' || icode == '\0') ? "" : std::string(1, icode);
  };

  // Convert gemmi altloc ('\\0' or ' ' = none) to coot convention ("")
  auto altloc_to_str = [](char alt) -> std::string {
    return (alt == '\0' || alt == ' ') ? "" : std::string(1, alt);
  };

  atom_spec_t a1(a1_addr.chain_name,
                 a1_addr.res_id.seqid.num.value,
                 icode_to_str(a1_addr.res_id.seqid.icode),
                 a1_addr.atom_name,
                 altloc_to_str(a1_addr.altloc));

  atom_spec_t a2(a2_addr.chain_name,
                 a2_addr.res_id.seqid.num.value,
                 icode_to_str(a2_addr.res_id.seqid.icode),
                 a2_addr.atom_name,
                 altloc_to_str(a2_addr.altloc));

  a1.model_number = model_number;
  a2.model_number = model_number;

  return {a1, a2};
}

} // namespace coot
