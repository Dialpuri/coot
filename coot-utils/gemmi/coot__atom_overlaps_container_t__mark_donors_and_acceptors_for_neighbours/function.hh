#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <unordered_map>
#include "geometry/hb-types.hh"
#include "coot-utils/atom-overlaps.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

void mark_donors_and_acceptors_for_neighbours_gemmi(
    const protein_geometry* geom_p,
    const std::vector<gemmi::Residue*>& neighbours,
    std::unordered_map<const gemmi::Atom*, hb_t>& atom_hb_types,
    bool& molecule_has_hydrogens
);

} // namespace coot