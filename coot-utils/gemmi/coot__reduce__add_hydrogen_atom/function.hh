#pragma once
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace reduce {

// Helper to trim whitespace from both ends of a string
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

/// Add a hydrogen atom to a residue. If an atom with the same name and
/// alternative location already exists, update its coordinates.
/// Returns the added or updated atom.
inline gemmi::Atom *
add_hydrogen_atom_gemmi(const std::string& atom_name, const clipper::Coord_orth& pos,
                        float bf, const std::string& altconf,
                        gemmi::Residue& residue) {
    gemmi::Atom new_H;
    new_H.name = atom_name;
    new_H.element = gemmi::Element("H");
    new_H.pos = gemmi::Position(pos.x(), pos.y(), pos.z());
    new_H.occ = 1.0f;
    new_H.b_iso = bf;
    if (!altconf.empty() && altconf.length() > 0) {
        new_H.altloc = altconf[0];
    }

    // Check if atom already exists with same name
    // Use trimmed comparison to handle leading/trailing spaces
    std::string trimmed_atom_name = trim(atom_name);
    for (auto& atom : residue.atoms) {
        if (trim(atom.name) == trimmed_atom_name) {
            if (altconf.empty()) {
                // Match any atom with this name when altconf is empty
                atom.pos = new_H.pos;
                atom.occ = new_H.occ;
                atom.b_iso = new_H.b_iso;
                atom.name = atom_name;
                return &atom;
            }
            if (atom.altloc == new_H.altloc) {
                atom.pos = new_H.pos;
                atom.occ = new_H.occ;
                atom.b_iso = new_H.b_iso;
                atom.name = atom_name;
                return &atom;
            }
        }
    }

    // Add new atom
    residue.atoms.push_back(new_H);
    return &residue.atoms.back();
}

} // namespace reduce
} // namespace coot