#pragma once

#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

// Free-function port of coot::protein_geometry::get_group(mmdb::Residue*) to gemmi
// Calls the existing get_group(std::string) member on the real protein_geometry object
inline std::string get_group_gemmi(
    coot::protein_geometry &pg,
    const gemmi::Residue &residue
) {
    std::string res_name = residue.name;
    try {
        std::string group = pg.get_group(res_name);
        return group;
    }
    catch(const std::runtime_error& e) {
        std::cout << "ERROR:: pg::get_group" << res_name << " " << e.what() << '\n';
        return std::string("");
    }
}

} // namespace coot
