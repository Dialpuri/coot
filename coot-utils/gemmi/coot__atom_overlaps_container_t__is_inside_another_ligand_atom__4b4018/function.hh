#pragma once
#include <map>
#include <vector>
#include <string>
#include <clipper/clipper.h>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

class atom_overlaps_container_t {
public:
    // Original MMDB types replaced with gemmi equivalents
    gemmi::Model *model;
    int imol_enc;
    bool have_dictionary;
    bool molecule_has_hydrogens;
    gemmi::Residue *res_central;
    std::vector<gemmi::Residue *> neighbours;
    int udd_h_bond_type_handle;
    int udd_residue_index_handle;
    double probe_radius;
    bool ignore_water_contacts_flag;

    // maps from atom index to list of (atom, radius) pairs
    std::map<int, std::vector<std::pair<gemmi::Atom *, double>>> ligand_atom_neighbour_map;

    atom_overlaps_container_t() 
        : model(nullptr), imol_enc(-1), have_dictionary(false), 
          molecule_has_hydrogens(false), res_central(nullptr),
          udd_h_bond_type_handle(-1), udd_residue_index_handle(-1),
          probe_radius(1.0), ignore_water_contacts_flag(false) {}

    // Ported function - checks if a point is inside another ligand atom
    bool is_inside_another_ligand_atom_gemmi(int idx,
                                              const clipper::Coord_orth &dot_pt) const {
        bool r = false;

        if (idx >= 0) {
            std::map<int, std::vector<std::pair<gemmi::Atom *, double>>>::const_iterator it;
            it = ligand_atom_neighbour_map.find(idx);

            if (it != ligand_atom_neighbour_map.end()) {
                const std::vector<std::pair<gemmi::Atom *, double>> &v = it->second;
                for (size_t i = 0; i < v.size(); i++) {
                    // Get coordinates from gemmi atom
                    clipper::Coord_orth pt(v[i].first->pos.x, v[i].first->pos.y, v[i].first->pos.z);
                    double dist_sqrd = (dot_pt - pt).lengthsq();

                    const double &radius_other = v[i].second;
                    if (dist_sqrd < radius_other * radius_other) {
                        r = true;
                        break;
                    }
                }
            }
        }
        return r;
    }
};

// Helper function to convert gemmi atom to clipper::Coord_orth
inline clipper::Coord_orth co_gemmi(const gemmi::Atom* at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

} // namespace coot