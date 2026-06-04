#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "clipper/core/xmap.h"
#include "clipper/core/coords.h"
#include "coot-utils/coot-map-utils.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

inline double sum_density_for_atoms_in_residue_gemmi(
    gemmi::Model& model,
    const std::string& cid,
    const std::vector<std::string>& atom_names,
    const clipper::Xmap<float>& xmap) {

    double v = 0.0;
    gemmi::Residue* residue_p = coot::cid_to_residue_gemmi(cid, model);
    if (residue_p) {
        for (const gemmi::Atom& at : residue_p->atoms) {
            // gemmi atom.name is unpadded (e.g. "CA"). The MMDB GetAtomName() returns
            // the raw 4-char PDB field (columns 13-16), where the element symbol is
            // right-justified starting at position 1: e.g. " CA ", " N  ", " CG1".
            std::string atom_name = std::string(1, ' ') + at.name;
            atom_name.resize(4, ' ');
            if (std::find(atom_names.begin(), atom_names.end(), atom_name) != atom_names.end()) {
                clipper::Coord_orth pos = coot::co_gemmi(&at);
                float d = coot::util::density_at_point(xmap, pos);
                v += static_cast<double>(d);
            }
        }
    }
    return v;
}

} // namespace coot