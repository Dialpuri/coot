#pragma once

#include <gemmi/model.hpp>
#include <gemmi/select.hpp>
#include <string>
#include <vector>

namespace coot {
namespace util {

std::vector<std::string>
chains_in_atom_selection_gemmi(gemmi::Model& model, const std::string& atom_selection) {
  std::vector<std::string> result;

  gemmi::Selection sel(atom_selection);
  auto filtered_chains = sel.chains(model);

  for (const auto& chain : filtered_chains) {
    result.push_back(chain.name);
  }

  return result;
}

} // namespace util
} // namespace coot