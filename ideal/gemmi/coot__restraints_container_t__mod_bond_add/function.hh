#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

namespace coot {

// Simplified restraints_container_t for gemmi port
struct restraints_container_t {
    std::vector<std::vector<int>> restraints_vec;
    std::vector<std::set<int>> bonded_atom_indices;

    // Helpers (no _gemmi port exists — translated inline)
    std::vector<bool> make_fixed_flags(int idx_1, int idx_2);
    void add(int type, int idx_1, int idx_2,
             const std::vector<bool>& fixed_flags,
             double dist, double esd, double dummy);

    // Gemmi port of mod_bond_add
    void mod_bond_add_gemmi(const struct chem_mod_bond& mod_bond,
                            gemmi::Residue& residue);
};

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
          new_value_dist_esd(new_value_dist_) {}
};

// BOND_RESTRAINT enum value (matches coot::simple_restraint.hh)
constexpr int BOND_RESTRAINT = 0;

// ---------- inline definitions ----------

inline std::vector<bool> restraints_container_t::make_fixed_flags(int, int) {
    return std::vector<bool>(2, false);
}

inline void restraints_container_t::add(int type, int idx_1, int idx_2,
                                         const std::vector<bool>&,
                                         double, double, double) {
    std::vector<int> restraint;
    restraint.push_back(type);
    restraint.push_back(idx_1);
    restraint.push_back(idx_2);
    restraints_vec.push_back(restraint);
}

// Trim leading and trailing whitespace from a PDB atom name
// (gemmi stores 4-char padded names; MMDB stores null-terminated names)
inline std::string trim_atom_name(std::string name) {
    auto lo = name.find_first_not_of(' ');
    if (lo == std::string::npos) return "";
    auto hi = name.find_last_not_of(' ');
    return name.substr(lo, hi - lo + 1);
}

inline void restraints_container_t::mod_bond_add_gemmi(
    const chem_mod_bond& mod_bond,
    gemmi::Residue& residue)
{
    int index_1 = -1, index_2 = -1;
    size_t n_atoms = residue.atoms.size();

    for (size_t iat_1 = 0; iat_1 < n_atoms; iat_1++) {
        std::string name_1 = trim_atom_name(residue.atoms[iat_1].name);
        if (name_1 == mod_bond.atom_id_1) {
            for (size_t iat_2 = 0; iat_2 < n_atoms; iat_2++) {
                std::string name_2 = trim_atom_name(residue.atoms[iat_2].name);
                if (name_2 == mod_bond.atom_id_2) {
                    char alt_1 = residue.atoms[iat_1].altloc;
                    char alt_2 = residue.atoms[iat_2].altloc;
                    // In gemmi, ' ' (space) means "no alt loc", equivalent to MMDB ""
                    if (alt_1 == ' ' || alt_2 == ' ' || alt_1 == alt_2) {
                        index_1 = static_cast<int>(iat_1);
                        index_2 = static_cast<int>(iat_2);

                        // Ensure bonded_atom_indices is large enough
                        int max_idx = std::max(index_1, index_2);
                        if (static_cast<size_t>(max_idx + 1) > bonded_atom_indices.size()) {
                            bonded_atom_indices.resize(static_cast<size_t>(max_idx + 1));
                        }

                        bonded_atom_indices[index_1].insert(index_2);
                        bonded_atom_indices[index_2].insert(index_1);

                        std::vector<bool> fixed_flags =
                            make_fixed_flags(index_1, index_2);

                        add(BOND_RESTRAINT, index_1, index_2,
                            fixed_flags,
                            mod_bond.new_value_dist,
                            mod_bond.new_value_dist_esd,
                            1.2);

                        // Duplicate inserts from original (no-op since set, but matches original)
                        bonded_atom_indices[index_1].insert(index_2);
                        bonded_atom_indices[index_2].insert(index_1);
                    }
                }
            }
        }
    }
}

} // namespace coot