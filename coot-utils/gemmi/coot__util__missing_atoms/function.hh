#pragma once

#include <string>
#include <vector>
#include <map>

#include <gemmi/model.hpp>

#include "coot-utils/coot-coord-extras.hh"
#include "geometry/protein-geometry.hh"

namespace coot {
namespace util {

struct missing_atom_info_gemmi {
    std::vector<std::string> residues_with_no_dictionary;
    std::vector<gemmi::Residue*> residues_with_missing_atoms;
    std::vector<std::pair<gemmi::Residue*, std::vector<gemmi::Atom*>>> atoms_in_coords_but_not_in_dict;
    std::map<gemmi::Residue*, std::vector<std::string>> residue_missing_atom_names_map;

    missing_atom_info_gemmi(
        const std::vector<std::string>& _residues_with_no_dictionary = {},
        const std::vector<gemmi::Residue*>& _residues_with_missing_atoms = {},
        const std::vector<std::pair<gemmi::Residue*, std::vector<gemmi::Atom*>>>& _atoms_in_coords_but_not_in_dict = {}
    )
        : residues_with_no_dictionary(_residues_with_no_dictionary),
          residues_with_missing_atoms(_residues_with_missing_atoms),
          atoms_in_coords_but_not_in_dict(_atoms_in_coords_but_not_in_dict) {}
};

missing_atom_info_gemmi missing_atoms_gemmi(
    gemmi::Structure* st,
    bool do_missing_hydrogen_atoms_flag,
    coot::protein_geometry* geom_p
);

} // namespace util
} // namespace coot
