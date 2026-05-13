#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace pli {

// Port of MMDB pli::is_a_metal to gemmi.
// Returns true if the residue name matches a known metal ion.
bool is_a_metal_gemmi(const gemmi::Residue& residue) {
  const std::string& res_name = residue.name;
  return res_name == "MG" || res_name == "CA" || res_name == "MN" ||
         res_name == "FE" || res_name == "K" || res_name == "NA" ||
         res_name == "CO" || res_name == "NI" || res_name == "CU" ||
         res_name == "ZN" || res_name == "RU" || res_name == "PT" ||
         res_name == "AU" || res_name == "AG";
}

} // namespace pli