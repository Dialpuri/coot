#pragma once
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>

namespace coot { namespace reduce {

// ── helpers (also used externally) ──────────────────────────

inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last  = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

inline std::vector<std::string> get_residue_alt_confs_gemmi(const gemmi::Residue& residue) {
    std::vector<std::string> v;
    for (const auto& atom : residue.atoms) {
        char alt = atom.altloc;
        std::string alt_str = (alt == ' ' || alt == '\0') ? "" : std::string(1, alt);
        bool found = false;
        for (const auto& s : v) { if (s == alt_str) { found = true; break; } }
        if (!found) v.push_back(alt_str);
    }
    return v;
}

inline const gemmi::Atom* find_atom_in_residue(const gemmi::Residue& res,
                                               const std::string& name,
                                               const std::string& altconf) {
    std::string trimmed_name = trim(name);
    for (const auto& atom : res.atoms) {
        if (trim(atom.name) != trimmed_name) continue;
        if (!altconf.empty()) {
            char alt = atom.altloc;
            std::string alt_str = (alt == ' ' || alt == '\0') ? "" : std::string(1, alt);
            if (alt_str != altconf) continue;
        }
        return &atom;
    }
    return nullptr;
}

inline gemmi::Atom* add_hydrogen_atom_gemmi(const std::string& atom_name,
                                            const clipper::Coord_orth& pos,
                                            float bf,
                                            const std::string& altconf,
                                            gemmi::Residue& residue) {
    gemmi::Atom new_H;
    new_H.name    = atom_name;
    new_H.element = gemmi::Element("H");
    new_H.pos     = gemmi::Position(pos.x(), pos.y(), pos.z());
    new_H.occ     = 1.0f;
    new_H.b_iso   = bf;
    if (!altconf.empty() && altconf.length() > 0) {
        new_H.altloc = altconf[0];
    }

    std::string trimmed_atom_name = trim(atom_name);
    for (auto& atom : residue.atoms) {
        if (trim(atom.name) != trimmed_atom_name) continue;
        if (altconf.empty()) {
            atom.pos   = new_H.pos;
            atom.occ   = new_H.occ;
            atom.b_iso = new_H.b_iso;
            atom.name  = atom_name;
            return &atom;
        }
        if (atom.altloc == new_H.altloc) {
            atom.pos   = new_H.pos;
            atom.occ   = new_H.occ;
            atom.b_iso = new_H.b_iso;
            atom.name  = atom_name;
            return &atom;
        }
    }
    residue.atoms.push_back(new_H);
    return &residue.atoms.back();
}

// ── add_xH_H (no _gemmi port exists — inline translation) ───

inline bool add_xH_H_gemmi(const std::string& H_at_name,
                           const std::string& first_neighb,
                           const std::string& second_neighb,
                           const std::string& third_neighb,
                           double bond_length,
                           double ang_deg,
                           double torsion_deg,
                           gemmi::Residue* residue_p) {
    if (!residue_p) return false;

    double ang_rad     = ang_deg     * M_PI / 180.0;
    double torsion_rad = torsion_deg * M_PI / 180.0;

    std::vector<std::string> alts = get_residue_alt_confs_gemmi(*residue_p);

    for (size_t i = 0; i < alts.size(); i++) {
        const gemmi::Atom* at_1 = find_atom_in_residue(*residue_p, first_neighb, alts[i]);
        if (!at_1) continue;
        const gemmi::Atom* at_2 = find_atom_in_residue(*residue_p, second_neighb, alts[i]);
        if (!at_2) continue;
        const gemmi::Atom* at_3 = find_atom_in_residue(*residue_p, third_neighb, alts[i]);
        if (!at_3) continue;

        clipper::Coord_orth p1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
        clipper::Coord_orth p2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
        clipper::Coord_orth p3(at_3->pos.x, at_3->pos.y, at_3->pos.z);

        clipper::Coord_orth new_pos(p1, p2, p3, bond_length, ang_rad, torsion_rad);

        float bf = 0;
        if (at_1->b_iso > 0) bf = at_1->b_iso;
        else if (at_2->b_iso > 0) bf = at_2->b_iso;
        else if (at_3->b_iso > 0) bf = at_3->b_iso;

        add_hydrogen_atom_gemmi(H_at_name, new_pos, bf, alts[i], *residue_p);
        return true;
    }
    return false;
}

// ── public add_OH_H ─────────────────────────────────────────

inline void add_OH_H_gemmi(const std::string &H_at_name,
                           const std::string &first_neighb,
                           const std::vector<std::string> &second_neighb_vec,
                           const std::map<std::string, std::vector<std::string>> &third_neighb_map,
                           double bond_length,
                           double ang_deg,
                           double torsion_deg,
                           gemmi::Residue *residue_p) {
    if (second_neighb_vec.size() > 0) {
        std::string second = second_neighb_vec[0];
        auto it = third_neighb_map.find(second);
        if (it != third_neighb_map.end()) {
            std::vector<std::string> thirds = it->second;
            if (thirds.size() > 0) {
                std::string third = thirds[0];
                add_xH_H_gemmi(H_at_name, first_neighb, second_neighb_vec[0], third,
                               bond_length, ang_deg, torsion_deg, residue_p);
            }
        } else {
            std::cout << "failed to find key " << second << " in thirds map" << std::endl;
        }
    }
}

}} // coot::reduce