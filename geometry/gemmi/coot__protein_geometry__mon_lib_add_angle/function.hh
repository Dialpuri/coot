#pragma once

#include <string>
#include <vector>
#include <utility>

// IMOL_ENC_AUTO value (from protein-geometry.hh)
namespace coot { namespace protein_geometry { constexpr int IMOL_ENC_AUTO = 0; } }

namespace coot {

// --- Minimal definitions mirroring coot types needed by mon_lib_add_angle ---

struct dict_angle_restraint_t {
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    double angle;
    double angle_esd;

    dict_angle_restraint_t(const std::string& a1, const std::string& a2, const std::string& a3,
                           double val, double esd)
        : atom_id_1(a1), atom_id_2(a2), atom_id_3(a3), angle(val), angle_esd(esd) {}
};

struct dict_chem_comp_t {
    std::string comp_id;
};

struct dictionary_residue_restraints_t {
    dict_chem_comp_t residue_info;
    std::vector<dict_angle_restraint_t> angle_restraint;

    dictionary_residue_restraints_t(const std::string& comp_id, int /*read_number*/)
        : residue_info(), angle_restraint()
    {
        residue_info.comp_id = comp_id;
    }
};

// --- gemmi port of mon_lib_add_angle ---

inline void mon_lib_add_angle_gemmi(
    const std::string& comp_id,
    int imol_enc,
    const std::string& atom_id_1,
    const std::string& atom_id_2,
    const std::string& atom_id_3,
    double value_angle,
    double value_angle_esd,
    std::vector<std::pair<int, dictionary_residue_restraints_t>>& dict_res_restraints,
    int& read_number)
{
    bool ifound = false;

    for (unsigned int i = 0; i < dict_res_restraints.size(); i++) {
        if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
            if (dict_res_restraints[i].first == imol_enc) {
                ifound = true;
                // Check for duplicate angle restraint
                bool already = false;
                for (unsigned int ia = 0; ia < dict_res_restraints[i].second.angle_restraint.size(); ia++) {
                    if (dict_res_restraints[i].second.angle_restraint[ia].atom_id_1 == atom_id_1 &&
                        dict_res_restraints[i].second.angle_restraint[ia].atom_id_2 == atom_id_2 &&
                        dict_res_restraints[i].second.angle_restraint[ia].atom_id_3 == atom_id_3) {
                        already = true;
                        break;
                    }
                }
                if (!already) {
                    dict_res_restraints[i].second.angle_restraint.push_back(
                        dict_angle_restraint_t(atom_id_1, atom_id_2, atom_id_3, value_angle, value_angle_esd));
                }
            }
        }
    }

    if (!ifound) {
        dictionary_residue_restraints_t rest(comp_id, read_number);
        std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
        dict_res_restraints.push_back(p);
        dict_res_restraints[dict_res_restraints.size() - 1].second.angle_restraint.push_back(
            dict_angle_restraint_t(atom_id_1, atom_id_2, atom_id_3, value_angle, value_angle_esd));
    }
}

} // namespace coot
