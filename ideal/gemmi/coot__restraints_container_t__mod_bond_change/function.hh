#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {

constexpr int BOND_RESTRAINT = 1;

struct simple_restraint {
    int restraint_type;
    int atom_index_1;
    int atom_index_2;
    float target_value;
    float sigma;
};

struct chem_mod_bond {
    std::string name;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string extra;
    float new_value_dist;
    float new_value_dist_esd;

    chem_mod_bond(std::string n, std::string a1, std::string a2, std::string e,
                  float d, float e_)
        : name(std::move(n)), atom_id_1(std::move(a1)), atom_id_2(std::move(a2)),
          extra(std::move(e)), new_value_dist(d), new_value_dist_esd(e_) {}
};

// Port of coot::restraints_container_t::mod_bond_change
//
// Finds bond restraints whose two atoms belong to target_residue and whose
// atom names match mod_bond.atom_id_1 / mod_bond.atom_id_2, then updates
// target_value and sigma.
void mod_bond_change_gemmi(std::vector<simple_restraint>& restraints_vec,
                           const std::vector<std::string>& atom_names,
                           const std::vector<gemmi::Residue*>& atom_residues,
                           const chem_mod_bond& mod_bond,
                           gemmi::Residue* target_residue);

}  // namespace coot