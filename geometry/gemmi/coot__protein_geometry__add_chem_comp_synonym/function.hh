#pragma once

#include <gemmi/cif.hpp>
#include <string>
#include <vector>

namespace coot {
namespace protein_geometry {

struct residue_name_synonym {
  std::string comp_id;
  std::string comp_alternative_id;
  std::string mod_id;

  residue_name_synonym() = default;
  residue_name_synonym(std::string c, std::string a, std::string m)
    : comp_id(std::move(c)), comp_alternative_id(std::move(a)), mod_id(std::move(m)) {}
};

inline void add_chem_comp_synonym_gemmi(
    const gemmi::cif::Loop& loop,
    std::vector<residue_name_synonym>& residue_name_synonyms) {

  int comp_id_idx = loop.find_tag("comp_id");
  int comp_alt_idx = loop.find_tag("comp_alternative_id");
  int mod_id_idx = loop.find_tag("mod_id");

  if (comp_id_idx == -1 || comp_alt_idx == -1 || mod_id_idx == -1)
    return;

  size_t n_rows = loop.length();
  for (size_t j = 0; j < n_rows; j++) {
    std::string comp_id = gemmi::cif::as_string(loop.val(j, comp_id_idx));
    std::string comp_alternative_id = gemmi::cif::as_string(loop.val(j, comp_alt_idx));
    std::string mod_id = gemmi::cif::as_string(loop.val(j, mod_id_idx));
    residue_name_synonyms.emplace_back(comp_id, comp_alternative_id, mod_id);
  }
}

} // namespace protein_geometry
} // namespace coot