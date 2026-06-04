#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

struct restraints_container_t {
   std::vector<std::pair<bool, gemmi::Residue*>> residues_vec;
   std::map<gemmi::Residue*, std::set<gemmi::Residue*>> fixed_neighbours_set;
   gemmi::Model* mol;
   std::vector<std::vector<int>> restraints_vec;

   void debug_sets_gemmi() const;
};

} // namespace coot