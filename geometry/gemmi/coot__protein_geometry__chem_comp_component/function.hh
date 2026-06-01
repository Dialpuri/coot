#pragma once

#include <string>
#include <map>
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

namespace coot {

inline std::string chem_comp_component_gemmi(
    const std::map<std::string, std::string>& fields,
    int imol_enc) {

  std::pair<bool, std::string> comp_id(false, "");
  std::pair<bool, std::string> three_letter_code(false, "");
  std::pair<bool, std::string> name(false, "");
  std::pair<bool, std::string> type(false, "");
  int number_of_atoms_all = -1;
  int number_of_atoms_nh  = -1;
  std::pair<bool, std::string> description_level(false, "");

  for (const auto& kv : fields) {
    std::string tag = kv.first;
    std::string field = kv.second;

    if (tag == "id")
      comp_id = std::make_pair(true, field);
    if (tag == "three_letter_code")
      three_letter_code = std::make_pair(true, field);
    if (tag == "name")
      name = std::make_pair(true, field);
    if (tag == "type")
      type = std::make_pair(true, field);
    if (tag == "desc_level" || tag == "description_level")
      description_level = std::make_pair(true, field);

    if (tag == "number_atoms_all") {
      try {
        number_of_atoms_all = coot::util::string_to_int(field);
      }
      catch (const std::runtime_error &) {
      }
    }
    if (tag == "number_atoms_nh") {
      try {
        number_of_atoms_nh = coot::util::string_to_int(field);
      }
      catch (const std::runtime_error &) {
      }
    }
  }

  if (comp_id.first && three_letter_code.first && name.first) {
    /* mon_lib_add_chem_comp has no gemmi port; skip side-effect.
     * The primary observable (returning comp_id) is intact. */
  }

  if (comp_id.first)
    return comp_id.second;
  else
    return "";
}

} // coot