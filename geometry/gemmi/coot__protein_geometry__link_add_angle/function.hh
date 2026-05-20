#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace protein_geometry {

// dict_link_angle_restraint_t equivalent
struct dict_link_angle_restraint_t {
    int atom_1_comp_id;
    int atom_2_comp_id;
    int atom_3_comp_id;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    double angle_val;
    double angle_esd_val;

    dict_link_angle_restraint_t(int a1, int a2, int a3,
                                const std::string& id1,
                                const std::string& id2,
                                const std::string& id3,
                                double val, double esd)
        : atom_1_comp_id(a1), atom_2_comp_id(a2), atom_3_comp_id(a3),
          atom_id_1(id1), atom_id_2(id2), atom_id_3(id3),
          angle_val(val), angle_esd_val(esd) {}

    static std::string pad_to_4c(const std::string& s) {
        if (s.size() >= 4) return s.substr(0, 4);
        // Pad on the right to make 4 characters (PDB-style atom names)
        // but add leading space if the name is short for right-alignment
        if (s.size() == 1) return "  " + s + " ";
        if (s.size() == 2) return " " + s + " ";
        if (s.size() == 3) return s + " ";
        return s + std::string(4 - s.size(), ' ');
    }
    std::string atom_id_1_4c() const { return pad_to_4c(atom_id_1); }
    std::string atom_id_2_4c() const { return pad_to_4c(atom_id_2); }
    std::string atom_id_3_4c() const { return pad_to_4c(atom_id_3); }
    double angle() const { return angle_val; }
    double angle_esd() const { return angle_esd_val; }
};

// dictionary_residue_link_restraints_t equivalent
struct dictionary_residue_link_restraints_t {
    std::string link_id;
    std::vector<dict_link_angle_restraint_t> link_angle_restraint;

    dictionary_residue_link_restraints_t(const std::string& lid)
        : link_id(lid), link_angle_restraint() {}
};

// Main function ported to gemmi
void link_add_angle_gemmi(const std::string& link_id,
                          int atom_1_comp_id,
                          int atom_2_comp_id,
                          int atom_3_comp_id,
                          const std::string& atom_id_1,
                          const std::string& atom_id_2,
                          const std::string& atom_id_3,
                          double value_dist,
                          double value_dist_esd,
                          std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints) {
    dict_link_angle_restraint_t lar(atom_1_comp_id,
                                    atom_2_comp_id,
                                    atom_3_comp_id,
                                    atom_id_1,
                                    atom_id_2,
                                    atom_id_3,
                                    value_dist,
                                    value_dist_esd);

    short int ifound = 0;
    for (size_t i = 0; i < dict_link_res_restraints.size(); i++) {
        if (dict_link_res_restraints[i].link_id == link_id) {
            ifound = 1;
            dict_link_res_restraints[i].link_angle_restraint.push_back(lar);
        }
    }

    if (!ifound) {
        dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
        dict_link_res_restraints[dict_link_res_restraints.size() - 1].link_angle_restraint.push_back(lar);
    }
}

} // namespace protein_geometry
} // namespace coot