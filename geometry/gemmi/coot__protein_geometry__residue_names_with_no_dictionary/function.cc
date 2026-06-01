#include "function.hh"

#include <set>
#include <algorithm>

namespace coot {

std::vector<std::string>
residue_names_with_no_dictionary_gemmi(
    const protein_geometry& geom,
    const gemmi::Structure* st,
    int imol_no) {

  std::vector<std::string> v;

  if (st && !st->models.empty()) {
    const gemmi::Model& model = st->models[0];
    std::set<std::string> already_tested_names;
    for (const gemmi::Chain& chain : model.chains) {
      for (const gemmi::Residue& residue : chain.residues) {
        std::string residue_name = residue.name;
        if (already_tested_names.find(residue_name) == already_tested_names.end()) {
          if (!geom.have_dictionary_for_residue_type_no_dynamic_add(residue_name, imol_no))
            if (std::find(v.begin(), v.end(), residue_name) == v.end())
              v.push_back(residue_name);
          already_tested_names.insert(residue_name);
        }
      }
    }
  }
  return v;
}

} // namespace coot