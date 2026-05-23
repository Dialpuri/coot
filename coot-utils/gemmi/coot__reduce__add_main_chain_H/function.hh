#pragma once
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

namespace coot {
namespace reduce {

// Helper function to trim whitespace
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

// Helper to get alternative conformations for a residue
inline std::vector<std::string> get_residue_alt_confs_gemmi(const gemmi::Residue& residue) {
    std::vector<std::string> v;
    for (const auto& atom : residue.atoms) {
        char alt = atom.altloc;
        std::string alt_str = (alt == ' ' || alt == '\0') ? "" : std::string(1, alt);
        
        bool ifound = false;
        for (size_t i = 0; i < v.size(); i++) {
            if (v[i] == alt_str) {
                ifound = true;
                break;
            }
        }
        if (!ifound)
            v.push_back(alt_str);
    }
    return v;
}

// Port of coot::reduce::add_hydrogen_atom
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

    std::string trimmed_atom_name = trim(atom_name);
    for (auto& atom : residue.atoms) {
        if (trim(atom.name) == trimmed_atom_name) {
            if (altconf.empty()) {
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

    residue.atoms.push_back(new_H);
    return &residue.atoms.back();
}

// Port of coot::reduce::add_main_chain_H from MMDB to gemmi
inline void
add_main_chain_H_gemmi(gemmi::Residue& residue, gemmi::Residue* residue_prev) {
    if (residue_prev) {
        std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);

        for (const std::string& alt_conf : alt_confs) {
            char alt_char = alt_conf.empty() ? ' ' : alt_conf[0];
            gemmi::Atom* at_ca = nullptr;
            gemmi::Atom* at_n = nullptr;
            gemmi::Atom* at_c_prev = nullptr;
            gemmi::Atom* at_o_prev = nullptr;

            for (auto& atom : residue.atoms) {
                std::string name_trimmed = trim(atom.name);
                if (name_trimmed == "CA" && (alt_char == ' ' || atom.altloc == alt_char)) {
                    at_ca = &atom;
                }
                if (name_trimmed == "N" && (alt_char == ' ' || atom.altloc == alt_char)) {
                    at_n = &atom;
                }
            }
            for (auto& atom : residue_prev->atoms) {
                std::string name_trimmed = trim(atom.name);
                if (name_trimmed == "C" && (alt_char == ' ' || atom.altloc == alt_char)) {
                    at_c_prev = &atom;
                }
                if (name_trimmed == "O" && (alt_char == ' ' || atom.altloc == alt_char)) {
                    at_o_prev = &atom;
                }
            }

            if (at_ca && at_n && at_c_prev && at_o_prev) {
                clipper::Coord_orth at_c_pos(at_c_prev->pos.x, at_c_prev->pos.y, at_c_prev->pos.z);
                clipper::Coord_orth at_o_pos(at_o_prev->pos.x, at_o_prev->pos.y, at_o_prev->pos.z);
                clipper::Coord_orth at_n_pos(at_n->pos.x, at_n->pos.y, at_n->pos.z);
                clipper::Coord_orth at_ca_pos(at_ca->pos.x, at_ca->pos.y, at_ca->pos.z);
                
                double bl = 0.86;
                double angle = clipper::Util::d2rad(125.0);
                clipper::Coord_orth H_pos(at_ca_pos, at_c_pos, at_n_pos, bl, angle, M_PI);
                
                float bf = at_n->b_iso;
                add_hydrogen_atom_gemmi(" H  ", H_pos, bf, alt_conf, residue);
            }
        }
    }
}

} // namespace reduce
} // namespace coot