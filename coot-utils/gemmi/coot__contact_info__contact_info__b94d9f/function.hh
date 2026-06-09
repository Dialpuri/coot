#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include "geometry/protein-geometry.hh"
#include "coot-utils/map-index.hh"

namespace coot {
namespace contact_info {

struct contacts_pair {
    int i;
    int j;
    contacts_pair(int i_in, int j_in) : i(i_in), j(j_in) {}
};

struct contact_info_result {
    std::vector<contacts_pair> contacts;
    int n_contacts() const { return static_cast<int>(contacts.size()); }
};

inline contact_info_result contact_info_gemmi(
    const std::vector<gemmi::Atom*>& atom_selection,
    const std::string& monomer_type,
    int imol,
    coot::protein_geometry* geom_p)
{
    contact_info_result result;

    std::pair<bool, coot::dictionary_residue_restraints_t> r =
        geom_p->get_monomer_restraints(monomer_type, imol);

    if (r.first) {
        std::map<std::string, coot::map_index_t> name_map;
        for (int i = 0; i < static_cast<int>(atom_selection.size()); i++) {
            std::string atom_name = coot::atom_id_mmdb_expand(atom_selection[i]->name);
            name_map[atom_name] = coot::map_index_t(i);
        }

        for (unsigned int ib = 0; ib < r.second.bond_restraint.size(); ib++) {
            std::string n_1 = r.second.bond_restraint[ib].atom_id_1_4c();
            std::string n_2 = r.second.bond_restraint[ib].atom_id_2_4c();
            coot::map_index_t ind_1 = name_map[n_1];
            coot::map_index_t ind_2 = name_map[n_2];
            if (ind_1.is_assigned() && ind_2.is_assigned()) {
                contacts_pair p(ind_1.index(), ind_2.index());
                result.contacts.push_back(p);
            }
        }
    }

    return result;
}

} // namespace contact_info
} // namespace coot
