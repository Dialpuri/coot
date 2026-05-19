#include "function.hh"
#include <iostream>

namespace coot {

atom_overlaps_container_t::atom_overlaps_container_t(
    const gemmi::Residue& res,
    const std::vector<const gemmi::Residue*>& neis,
    const gemmi::Model& mod,
    const void* geom)
    : central_residue(&res),
      neighbours(neis),
      geom_p(geom) {}

double atom_overlaps_container_t::get_vdw_radius_ligand_atom_gemmi(const gemmi::Atom* at) const {
    double r = 2.5;

    auto it = central_residue_atoms_vdw_radius_map.find(at);
    if (it == central_residue_atoms_vdw_radius_map.end()) {
        // we need to add it then

        // What's the energy type of Atom at?
        std::string te = type_energy(at);
        if (!te.empty()) {
            auto it_type = type_to_vdw_radius_map.find(te);
            if (it_type == type_to_vdw_radius_map.end()) {
                // didn't find it. so look it up and add it.
                if (geom_p)
                    r = type_energy_to_radius(te);
                else
                    r = type_energy_to_radius(te);  // Also compute when geom_p is null
                type_to_vdw_radius_map[te] = r;
            } else {
                r = it_type->second;
            }
            central_residue_atoms_vdw_radius_map[at] = r;
        } else {
            std::cout << "failed to find type-energy for atom " << at->name << std::endl;
        }
    } else {
        r = it->second;
    }

    return r;
}

std::string atom_overlaps_container_t::type_energy(const gemmi::Atom* at) const {
    // Simple implementation - in real coot this would use dictionary
    // For now, return the atom name as a placeholder
    return at->name;
}

double atom_overlaps_container_t::type_energy_to_radius(const std::string& te) const {
    // Simple implementation - in real coot this would use geometry information
    // Default radius for unknown types
    if (te == "O") return 1.52;
    if (te == "O1" || te == "O2" || te == "O3" || te == "O4" || te == "O5" || te == "O6" || te == "O7" || te == "O8" || te == "O9") return 1.52;
    return 2.5;
}

} // namespace coot