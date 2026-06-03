#pragma once

#include <map>
#include <string>
#include <iostream>
#include "geometry/protein-geometry.hh"
#include "geometry/gphl-chem-comp-info.hh"

namespace coot {

// Port of protein_geometry::gphl_chem_comp_info (private member function)
// Uses std::map<std::string,std::string> instead of mmdb::mmcif::PStruct
// Returns the built gphl_chem_comp_info_t (empty if no comp_id or unknown residue)
inline gphl_chem_comp_info_t gphl_chem_comp_info_gemmi(
    protein_geometry& geom,
    const std::map<std::string, std::string>& structure_data,
    int imol_enc)
{
    gphl_chem_comp_info_t gphl_chem_comp_info;

    for (const auto& [tag, field] : structure_data) {
        gphl_chem_comp_info.add(tag, field);
    }

    int idx = gphl_chem_comp_info.get_index("comp_id");
    if (idx >= 0) {
        const std::string& comp_id = gphl_chem_comp_info[idx].second;
        int idx_rest = geom.get_monomer_restraints_index(comp_id, imol_enc, true);
        if (idx_rest >= 0) {
            std::cout << "debug:: adding a gphl info for " << comp_id
                      << " of size " << gphl_chem_comp_info.info.size() << std::endl;
            return gphl_chem_comp_info;
        }
    }
    return gphl_chem_comp_info_t{};
}

} // namespace coot