#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <utility>
#include <vector>
#include <string>

namespace coot {
namespace db_strands {

static bool atom_name_matches(const std::string& name, const std::string& target) {
    // Trim whitespace and compare
    size_t start = name.find_first_not_of(" \t");
    if (start == std::string::npos)
        return false;
    size_t end = name.find_last_not_of(" \t");
    return name.substr(start, end - start + 1) == target;
}

std::pair<bool, clipper::RTop_orth>
orient_strand_on_z_gemmi(const std::vector<const gemmi::Residue*>& residues) {

    clipper::Mat33<double> m_dum(1, 0, 0, 0, 1, 0, 0, 0, 1);
    clipper::Coord_orth pt_dum(0, 0, 0);
    clipper::RTop_orth rtop(m_dum, pt_dum);
    bool stat = false;

    int n_residues = static_cast<int>(residues.size());

    // z_control_points logic inlined (no _gemmi port exists)
    // Each residue gets 3 z-points (one for N, CA, C)
    const double spacing = 3.8;
    std::vector<clipper::Coord_orth> z_points;
    if (n_residues > 0) {
        double z_start = -spacing * (n_residues - 1) / 2.0;
        for (int i = 0; i < n_residues; i++) {
            clipper::Coord_orth zp(0, 0, z_start + spacing * i);
            z_points.push_back(zp);
            z_points.push_back(zp);
            z_points.push_back(zp);
        }
    }

    std::vector<clipper::Coord_orth> atom_vec;

    for (const gemmi::Residue* res : residues) {
        for (const gemmi::Atom& atom : res->atoms) {
            std::string atom_name = atom.name;
            if (atom_name_matches(atom_name, "N")) {
                atom_vec.push_back(clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
            }
            if (atom_name_matches(atom_name, "CA")) {
                atom_vec.push_back(clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
            }
            if (atom_name_matches(atom_name, "C")) {
                atom_vec.push_back(clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
            }
        }
    }

    if (n_residues == 0) {
        // Empty selection — return success with identity transform (original behavior)
        stat = true;
    } else if (static_cast<int>(atom_vec.size()) != 3 * n_residues) {
        // Missing backbone atoms — fail
        stat = false;
    } else {
        rtop = clipper::RTop_orth(atom_vec, z_points);
        stat = true;
    }

    return {stat, rtop};
}

} // namespace db_strands
} // namespace coot