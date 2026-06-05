#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <algorithm>

namespace coot {

// chem_mod_bond — chemical modification bond descriptor
struct chem_mod_bond {
    std::string type;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string chain_id;
    double new_value_dist;
    double new_value_dist_esd;

    chem_mod_bond(std::string type_,
                  std::string atom_id_1_,
                  std::string atom_id_2_,
                  std::string chain_id_,
                  double new_value_dist_,
                  double new_value_dist_esd_)
        : type(std::move(type_)),
          atom_id_1(std::move(atom_id_1_)),
          atom_id_2(std::move(atom_id_2_)),
          chain_id(std::move(chain_id_)),
          new_value_dist(new_value_dist_),
          new_value_dist_esd(new_value_dist_esd_) {}
};

// Bond restraint type identifier (value 0 matches coot::BOND_RESTRAINT)
static constexpr int BOND_TYPE_ID = 0;

// Simplified atom_spec for gemmi port — carries residue pointer and name
struct gemmi_atom_spec {
    const gemmi::Residue *residue;
    std::string name;  // 4-char padded atom name (PDB format)
};

// Simplified restraints_container_t for gemmi port
struct restraints_container_t {
    std::vector<std::vector<int>> restraints_vec;  // each: {type, idx_1, idx_2}
    std::vector<gemmi_atom_spec> atoms;

    // Gemmi port of mod_bond_delete
    void mod_bond_delete_gemmi(const chem_mod_bond &mod_bond,
                               const gemmi::Residue *residue_p);
};

// Helper: trim whitespace from both ends of a string
static inline std::string trim_atom_name(const std::string &s) {
    auto start = std::find_if_not(s.begin(), s.end(), ::isspace);
    auto end   = std::find_if_not(s.rbegin(), s.rend(), ::isspace).base();
    if (start >= end) return "";
    return std::string(start, end);
}

// ---------- inline definition ----------

inline void restraints_container_t::mod_bond_delete_gemmi(
    const chem_mod_bond &mod_bond,
    const gemmi::Residue *residue_p)
{
    // Trim the mod_bond atom IDs to get bare names for comparison
    std::string target_1 = trim_atom_name(mod_bond.atom_id_1);
    std::string target_2 = trim_atom_name(mod_bond.atom_id_2);

    for (auto it = restraints_vec.begin(); it != restraints_vec.end(); ) {
        if (it->size() >= 3 && (*it)[0] == BOND_TYPE_ID) {
            int idx_1 = (*it)[1];
            int idx_2 = (*it)[2];
            if (idx_1 >= 0 && idx_1 < static_cast<int>(atoms.size()) &&
                idx_2 >= 0 && idx_2 < static_cast<int>(atoms.size())) {
                if (atoms[idx_1].residue == residue_p) {
                    if (atoms[idx_2].residue == residue_p) {
                        std::string name_1 = trim_atom_name(atoms[idx_1].name);
                        std::string name_2 = trim_atom_name(atoms[idx_2].name);
                        if (name_1 == target_1) {
                            if (name_2 == target_2) {
                                it = restraints_vec.erase(it);
                                continue;
                            }
                        }
                    }
                }
            }
        }
        ++it;
    }
}

} // namespace coot