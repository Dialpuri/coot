#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>
#include <vector>

namespace coot {
namespace util {

// gemmi port of create_mmdbmanager_from_points.
// Returns a gemmi::Structure (value type) instead of mmdb::Manager*.
inline gemmi::Structure create_mmdbmanager_from_points_gemmi(
    const std::vector<clipper::Coord_orth>& pts, float b_factor) {

    gemmi::Structure st;
    gemmi::Model& model = st.models.emplace_back();

    gemmi::Chain& chain = model.chains.emplace_back();
    chain.name = "A";

    for (unsigned int i = 0; i < pts.size(); i++) {
        gemmi::Residue res;
        res.name = "ALA";
        res.seqid.num = i;

        gemmi::Atom atom;
        atom.name = " CA ";
        atom.element = gemmi::Element("C");
        atom.pos = gemmi::Position(pts[i].x(), pts[i].y(), pts[i].z());
        atom.occ = 1.0f;
        atom.b_iso = b_factor;

        res.atoms.push_back(atom);
        chain.residues.push_back(res);
    }

    return st;
}

} // namespace util
} // namespace coot