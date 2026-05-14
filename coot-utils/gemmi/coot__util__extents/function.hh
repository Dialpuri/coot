#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <limits>

namespace coot { namespace util {

// Port of coot::util::extents from MMDB to gemmi
// Returns pair of min/max coordinates (clipper::Coord_orth) for all atoms in the structure
inline std::pair<clipper::Coord_orth, clipper::Coord_orth>
extents_gemmi(const gemmi::Structure& st) {
    if (st.models.empty()) {
        return std::make_pair(
            clipper::Coord_orth(0.0, 0.0, 0.0),
            clipper::Coord_orth(0.0, 0.0, 0.0)
        );
    }

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                for (const gemmi::Atom& atom : residue.atoms) {
                    double x = atom.pos.x;
                    double y = atom.pos.y;
                    double z = atom.pos.z;
                    if (x < min_x) min_x = x;
                    if (y < min_y) min_y = y;
                    if (z < min_z) min_z = z;
                    if (x > max_x) max_x = x;
                    if (y > max_y) max_y = y;
                    if (z > max_z) max_z = z;
                }
            }
        }
    }

    // If no atoms found, return zeros
    if (std::isfinite(min_x) && std::isfinite(min_y) && std::isfinite(min_z)) {
        return std::make_pair(
            clipper::Coord_orth(min_x, min_y, min_z),
            clipper::Coord_orth(max_x, max_y, max_z)
        );
    } else {
        return std::make_pair(
            clipper::Coord_orth(0.0, 0.0, 0.0),
            clipper::Coord_orth(0.0, 0.0, 0.0)
        );
    }
}

}} // namespace coot::util