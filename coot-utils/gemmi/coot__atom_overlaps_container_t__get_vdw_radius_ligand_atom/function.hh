#pragma once
#include <string>
#include <map>
#include <gemmi/model.hpp>

namespace coot {

struct atom_overlaps_container_t {
    // Member variables that need to be present
    const gemmi::Residue* central_residue;
    std::vector<const gemmi::Residue*> neighbours;
    const void* geom_p;  // placeholder for protein_geometry pointer
    mutable std::map<const gemmi::Atom*, double> central_residue_atoms_vdw_radius_map;
    mutable std::map<std::string, double> type_to_vdw_radius_map;
    
    // Constructor - matching test signature
    atom_overlaps_container_t(const gemmi::Residue& res,
                              const std::vector<const gemmi::Residue*>& neis,
                              const gemmi::Model& mod,
                              const void* geom = nullptr);
    
    // The ported function - takes gemmi::Atom pointer
    double get_vdw_radius_ligand_atom_gemmi(const gemmi::Atom* at) const;
    
    // Helper functions
    std::string type_energy(const gemmi::Atom* at) const;
    double type_energy_to_radius(const std::string& te) const;
};

} // namespace coot