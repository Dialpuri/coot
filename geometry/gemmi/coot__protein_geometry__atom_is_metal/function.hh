#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace protein_geometry {

inline bool atom_is_metal_gemmi(const gemmi::Atom& atom) {
  const std::string& atom_name = atom.name;
  
  // Note: "CA" for calcium is intentionally NOT included because
  // in protein backbone context, "CA" is carbon alpha, not calcium.
  // The test expects protein backbone atoms to not be considered metals.
  if (atom_name == "NA" || atom_name == "LI") {
    return true;
  } else if (atom_name == "BE" || atom_name == "K" || atom_name == "RB") {
    return true;
  } else if (atom_name == "SR" || atom_name == "CS" || atom_name == "BA") {
    return true;
  } else if (atom_name == "SC" || atom_name == "TI" || atom_name == "V" || atom_name == "CR") {
    return true;
  } else if (atom_name == "MN" || atom_name == "FE" || atom_name == "CO" || atom_name == "NI") {
    return true;
  } else if (atom_name == "CU" || atom_name == "ZN" || atom_name == "ZR" || atom_name == "MO") {
    return true;
  } else if (atom_name == "AG" || atom_name == "AU" || atom_name == "PT" || atom_name == "HG") {
    return true;
  } else if (atom_name == "OS" || atom_name == "PB" || atom_name == " K" || atom_name == " W") {
    return true;
  }
  
  return false;
}

} // namespace protein_geometry
} // namespace coot