#pragma once

#include <string>
#include <vector>
#include <gemmi/cif.hpp>

namespace coot {

struct dictionary_residue_restraints_t {
  void write_cif_pdbx_chem_comp_descriptor_gemmi(gemmi::cif::Document* doc) const {
    // Placeholder implementation — matches original empty body
  }
};

} // namespace coot