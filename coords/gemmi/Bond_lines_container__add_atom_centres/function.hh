#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <map>
#include <string>

namespace coot {

// Gemmi port of Bond_lines_container::add_atom_centres
// Signature adapted to use gemmi::Structure instead of MMDB types
void add_atom_centres_gemmi(
    std::vector<gemmi::Vec3>& atom_centres,
    std::vector<int>& atom_centres_colour,
    const gemmi::Structure& st,
    int atom_colour_type,
    int model_number,
    bool do_bonds_to_hydrogens,
    int n_selected_atoms_hint,
    std::map<std::string, int>* atom_colour_map_p = nullptr);

}  // namespace coot