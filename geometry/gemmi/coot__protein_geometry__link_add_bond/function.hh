#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct dict_link_bond_restraint_t {
    int atom_1_comp_id;
    int atom_2_comp_id;
    std::string atom_id_1_;
    std::string atom_id_2_;
    double value_dist;
    double value_dist_esd;

    dict_link_bond_restraint_t(int a1, int a2,
                               const std::string& id1,
                               const std::string& id2,
                               double val, double esd)
        : atom_1_comp_id(a1), atom_2_comp_id(a2),
          atom_id_1_(id1), atom_id_2_(id2),
          value_dist(val), value_dist_esd(esd) {}

    std::string atom_id_1() const { return atom_id_1_; }
    std::string atom_id_2() const { return atom_id_2_; }
    double dist() const { return value_dist; }
    double esd() const { return value_dist_esd; }
};

struct dictionary_residue_link_restraints_t {
    std::string link_id;
    std::vector<dict_link_bond_restraint_t> link_bond_restraint;

    dictionary_residue_link_restraints_t(const std::string& lid)
        : link_id(lid), link_bond_restraint() {}
    dictionary_residue_link_restraints_t() : link_id("") {}
    bool empty() const { return link_id.empty(); }
};

void link_add_bond_gemmi(const std::string& link_id,
                         int atom_1_comp_id,
                         int atom_2_comp_id,
                         const std::string& atom_id_1,
                         const std::string& atom_id_2,
                         double value_dist,
                         double value_dist_esd,
                         std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints) {
    dict_link_bond_restraint_t lbr(atom_1_comp_id,
                                   atom_2_comp_id,
                                   atom_id_1,
                                   atom_id_2,
                                   value_dist,
                                   value_dist_esd);

    short int ifound = 0;
    for (size_t i = 0; i < dict_link_res_restraints.size(); i++) {
        if (dict_link_res_restraints[i].link_id == link_id) {
            ifound = 1;
            dict_link_res_restraints[i].link_bond_restraint.push_back(lbr);
        }
    }

    if (!ifound) {
        dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
        dict_link_res_restraints[dict_link_res_restraints.size() - 1].link_bond_restraint.push_back(lbr);
    }
}

} // namespace coot