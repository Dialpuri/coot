#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

namespace coot {

class protein_geometry;
class dictionary_residue_restraints_t;

double get_vdw_radius_neighb_atom_gemmi(
    const gemmi::Atom* at,
    unsigned int idx_res,
    const std::vector<dictionary_residue_restraints_t>& neighb_dictionaries,
    const protein_geometry* geom_p,
    std::map<const gemmi::Atom*, double>& neighbour_atoms_vdw_radius_map,
    std::map<std::string, double>& type_to_vdw_radius_map);

}  // namespace coot
