#pragma once

#include <gemmi/cifdoc.hpp>
#include <string>

namespace coot {
namespace protein_geometry {

/// Gemmi port of pdbx_chem_comp_model: extract the chemical component ID
/// from an mmCIF block (the _chem_comp.id tag).
inline std::string
pdbx_chem_comp_model_gemmi(const gemmi::cif::Block& block, int /*imol_enc*/) {
  const std::string* val = block.find_value("_chem_comp.id");
  if (val)
    return *val;
  return "";
}

} // namespace protein_geometry
} // namespace coot