#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

namespace coot {

// ----------------------------------------------------------------
// chem_mod_plane — chemical modification plane descriptor (gemmi port)
// ----------------------------------------------------------------
struct chem_mod_plane_gemmi {
    std::string name;
    std::string action;
    std::vector<std::pair<std::string, double>> atom_id_esd;

    chem_mod_plane_gemmi(std::string name_, std::string action_)
        : name(std::move(name_)), action(std::move(action_)) {}

    void add_atom_to_plane(std::string atom_id, double esd) {
        atom_id_esd.push_back(std::make_pair(std::move(atom_id), esd));
    }
};

// ----------------------------------------------------------------
// PLANE_RESTRAINT — restraint type constant (matches simple_restraint.hh)
// ----------------------------------------------------------------
constexpr int PLANE_RESTRAINT = 4;

// ----------------------------------------------------------------
// restraints_container_t_gemmi — self-contained gemmi version
// ----------------------------------------------------------------
struct restraints_container_t_gemmi {
    std::vector<std::vector<int>> restraints_vec;

    int size() const { return static_cast<int>(restraints_vec.size()); }

    // ----------------------------------------------------------------
    // make_fixed_flags — returns n booleans (one per atom), all false
    // ----------------------------------------------------------------
    std::vector<bool> make_fixed_flags(const std::vector<int> &indices) const {
        return std::vector<bool>(indices.size(), false);
    }

    // ----------------------------------------------------------------
    // add_plane — adds a plane restraint (type=4) if not duplicate
    // ----------------------------------------------------------------
    void add_plane(const std::vector<std::pair<int, double>> &position_sigma_indices,
                   const std::vector<bool> &fixed_flags) {
        // Check for duplicate
        for (const auto &r : restraints_vec) {
            if (static_cast<int>(r.size()) ==
                    static_cast<int>(position_sigma_indices.size() + 1) &&
                r[0] == PLANE_RESTRAINT) {
                bool match = true;
                for (size_t k = 0; k < position_sigma_indices.size(); ++k) {
                    if (r[k + 1] != position_sigma_indices[k].first) {
                        match = false;
                        break;
                    }
                }
                if (match) return; // duplicate — skip
            }
        }

        // Store restraint: [type, idx1, idx2, idx3, idx4, ...]
        std::vector<int> restraint;
        restraint.push_back(PLANE_RESTRAINT);
        for (const auto &ps : position_sigma_indices) {
            restraint.push_back(ps.first);
        }
        restraints_vec.push_back(restraint);
    }

    // ----------------------------------------------------------------
    // mod_plane_add_gemmi — the ported function
    // ----------------------------------------------------------------
    void mod_plane_add_gemmi(const struct chem_mod_plane_gemmi &mod_plane,
                             gemmi::Residue *residue_p);
};

// ----------------------------------------------------------------
// Helper: trim PDB atom name (removes padding spaces)
// ----------------------------------------------------------------
inline std::string trim_atom_name(std::string name) {
    auto lo = name.find_first_not_of(' ');
    if (lo == std::string::npos) return "";
    auto hi = name.find_last_not_of(' ');
    return name.substr(lo, hi - lo + 1);
}

// ----------------------------------------------------------------
// mod_plane_add_gemmi implementation
// ----------------------------------------------------------------
inline void restraints_container_t_gemmi::mod_plane_add_gemmi(
    const chem_mod_plane_gemmi &mod_plane,
    gemmi::Residue *residue_p)
{
    std::map<std::string, std::vector<int>> pos; // we worry about alt confs.

    for (unsigned int i = 0; i < mod_plane.atom_id_esd.size(); i++) {
        for (int iat = 0; iat < int(residue_p->atoms.size()); iat++) {
            std::string atom_name =
                trim_atom_name(std::string(residue_p->atoms[iat].name));
            if (atom_name == mod_plane.atom_id_esd[i].first) {
                // Use local residue atom index as the atom index
                int atom_index = iat;
                std::string altconf(1, residue_p->atoms[iat].altloc);
                pos[altconf].push_back(atom_index);
            }
        }
    }

    // iterate through all the alt confs (almost certainly only one)
    for (auto it = pos.begin(); it != pos.end(); it++) {
        const std::vector<int> &position_indices = it->second;

        if (position_indices.size() > 3) {
            std::vector<std::pair<int, double>> position_sigma_indices;
            for (unsigned int ii = 0; ii < position_indices.size(); ii++)
                position_sigma_indices.push_back(
                    std::pair<int, double>(position_indices[ii], 0.02));

            std::vector<bool> fixed_flags =
                make_fixed_flags(position_indices);
            add_plane(position_sigma_indices, fixed_flags);
        }
    }
}

} // namespace coot
