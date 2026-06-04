#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {
namespace util {

bool is_nucleotide_by_dict_dynamic_add_gemmi(
    gemmi::Residue* residue_p,
    coot::protein_geometry* geom_p
);

} // namespace util
} // namespace coot