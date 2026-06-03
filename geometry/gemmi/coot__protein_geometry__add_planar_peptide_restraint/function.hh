#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace protein_geometry {

// dict_link_plane_restraint_t equivalent
struct dict_link_plane_restraint_t {
    std::string plane_id;
    double dist_esd;
    std::vector<std::string> atom_ids;
    std::vector<int> atom_comp_ids;

    dict_link_plane_restraint_t(const std::string& pid, double desd)
        : plane_id(pid), dist_esd(desd) {}

    int n_atoms() const { return static_cast<int>(atom_ids.size()); }
    void set_dist_esd(double desd) { dist_esd = desd; }
};

// dictionary_residue_link_restraints_t equivalent
struct dictionary_residue_link_restraints_t {
    std::string link_id;
    std::vector<dict_link_plane_restraint_t> link_plane_restraint;

    dictionary_residue_link_restraints_t(const std::string& lid)
        : link_id(lid), link_plane_restraint() {}
};

// Replicate the private link_add_plane logic as a free function
void link_add_plane_gemmi(const std::string& link_id,
                          const std::string& atom_id,
                          const std::string& plane_id,
                          int atom_comp_id,
                          double dist_esd,
                          std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints) {
    dict_link_plane_restraint_t plr(plane_id, dist_esd);
    plr.atom_ids.push_back(atom_id);
    plr.atom_comp_ids.push_back(atom_comp_id);

    short int ifound = 0;
    for (size_t i = 0; i < dict_link_res_restraints.size(); i++) {
        if (dict_link_res_restraints[i].link_id == link_id) {
            // Check if this plane_id already exists in this link
            bool plane_found = false;
            for (size_t ip = 0; ip < dict_link_res_restraints[i].link_plane_restraint.size(); ip++) {
                if (dict_link_res_restraints[i].link_plane_restraint[ip].plane_id == plane_id) {
                    plane_found = true;
                    // Check if this atom is already in the plane restraint
                    int found_atom_index = -1;
                    for (int i_rest_at = 0; i_rest_at < dict_link_res_restraints[i].link_plane_restraint[ip].n_atoms(); i_rest_at++) {
                        if (dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids[i_rest_at] == atom_id) {
                            if (dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids[i_rest_at] == atom_comp_id) {
                                found_atom_index = i_rest_at;
                                break;
                            }
                        }
                    }
                    if (found_atom_index != -1) {
                        // Replace the atom
                        dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids[found_atom_index] = atom_id;
                        dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids[found_atom_index] = atom_comp_id;
                    } else {
                        // Append the atom
                        dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids.push_back(atom_id);
                        dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids.push_back(atom_comp_id);
                    }
                    ifound = 1;
                    break;
                }
            }
            if (!plane_found) {
                dict_link_res_restraints[i].link_plane_restraint.push_back(plr);
            }
            ifound = 1;
        }
    }

    if (!ifound) {
        dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
        dict_link_res_restraints[dict_link_res_restraints.size() - 1].link_plane_restraint.push_back(plr);
    }
}

// Port of add_planar_peptide_restraint to gemmi (free function)
void add_planar_peptide_restraint_gemmi(std::vector<dictionary_residue_link_restraints_t>& dict_link_res_restraints,
                                        bool& planar_peptide_restraint_state) {
    std::string plane_id = "plane-5-atoms";
    double dist_esd = 0.08;

    std::vector<std::pair<int, std::string>> v;
    v.push_back(std::pair<int, std::string>(1, "CA"));
    v.push_back(std::pair<int, std::string>(1, "C"));
    v.push_back(std::pair<int, std::string>(1, "O"));
    v.push_back(std::pair<int, std::string>(2, "N"));
    v.push_back(std::pair<int, std::string>(2, "CA"));

    for (unsigned int i = 0; i < v.size(); i++) {
        link_add_plane_gemmi("TRANS",  v[i].second, plane_id, v[i].first, dist_esd, dict_link_res_restraints);
        link_add_plane_gemmi("PTRANS", v[i].second, plane_id, v[i].first, dist_esd, dict_link_res_restraints);
    }
    planar_peptide_restraint_state = true;
}

} // namespace protein_geometry
} // namespace coot