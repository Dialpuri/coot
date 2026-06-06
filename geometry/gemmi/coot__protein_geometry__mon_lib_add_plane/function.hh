#pragma once

#include <string>
#include <vector>
#include <utility>

// IMOL_ENC_AUTO value (from protein-geometry.hh)
namespace coot { namespace protein_geometry { constexpr int IMOL_ENC_AUTO = 0; } }

namespace coot {

// --- Minimal definitions mirroring coot types needed by mon_lib_add_plane ---

struct dict_plane_restraint_t {
    std::string plane_id;
    double dist_esd;
    std::vector<std::string> atom_ids;

    dict_plane_restraint_t(const std::string& pid, const std::string& atom_id, double desd)
        : plane_id(pid), dist_esd(desd)
    {
        atom_ids.push_back(atom_id);
    }

    void push_back_atom(const std::string& atom_id, double /*dist_esd*/) {
        atom_ids.push_back(atom_id);
    }
};

struct dict_chem_comp_t {
    std::string comp_id;
};

struct dictionary_residue_restraints_t {
    dict_chem_comp_t residue_info;
    std::vector<dict_plane_restraint_t> plane_restraint;

    dictionary_residue_restraints_t(const std::string& comp_id, int /*read_number*/)
        : residue_info(), plane_restraint()
    {
        residue_info.comp_id = comp_id;
    }
};

// --- gemmi port of mon_lib_add_plane ---

inline void mon_lib_add_plane_gemmi(
    const std::string& comp_id,
    int imol_enc,
    const std::string& plane_id,
    const std::string& atom_id,
    double dist_esd,
    std::vector<std::pair<int, dictionary_residue_restraints_t>>& dict_res_restraints,
    int& read_number)
{
    bool ifound = false;

    for (unsigned int i = 0; i < dict_res_restraints.size(); i++) {
        if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
            if (dict_res_restraints[i].first == imol_enc) {
                for (unsigned int ip = 0; ip < dict_res_restraints[i].second.plane_restraint.size(); ip++) {
                    if (dict_res_restraints[i].second.plane_restraint[ip].plane_id == plane_id) {
                        ifound = true;
                        dict_res_restraints[i].second.plane_restraint[ip].push_back_atom(atom_id, dist_esd);
                        break;
                    }
                }
                if (!ifound) {
                    dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
                    dict_res_restraints[i].second.plane_restraint.push_back(res);
                    ifound = true;
                }
            }
        }
    }

    if (!ifound) {
        dictionary_residue_restraints_t rest(comp_id, read_number);
        std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
        dict_res_restraints.push_back(p);
        dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
        dict_res_restraints[dict_res_restraints.size() - 1].second.plane_restraint.push_back(res);
    }
}

} // namespace coot