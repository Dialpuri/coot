#pragma once

#include <string>
#include <gemmi/model.hpp>

// Forward declarations
namespace coot {
class protein_geometry;
struct atom_name_quad;
}

namespace coot {

int rotate_around_bond_gemmi(
    const std::string& residue_cid,
    const std::string& alt_conf,
    coot::atom_name_quad quad,
    double torsion_angle,
    coot::protein_geometry& geom,
    gemmi::Model& model);

} // namespace coot