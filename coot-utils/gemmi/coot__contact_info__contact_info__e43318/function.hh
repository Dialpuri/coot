#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>

#include "geometry/protein-geometry.hh"
#include "coot-utils/map-index.hh"

namespace coot {
namespace contact_info {

struct contacts_pair {
    int i;
    int j;
    contacts_pair(int i_in, int j_in) : i(i_in), j(j_in) {}
};

struct contact_info_result_gemmi {
    std::vector<contacts_pair> contacts;
    int n_contacts() const { return static_cast<int>(contacts.size()); }
};

inline contact_info_result_gemmi contact_info_gemmi(
    const std::vector<gemmi::Atom*>& atom_selection,
    const std::string& alt_conf,
    const std::string& monomer_type,
    int imol,
    coot::protein_geometry* geom_p)
{
    contact_info_result_gemmi result;

    std::pair<bool, coot::dictionary_residue_restraints_t> r =
        geom_p->get_monomer_restraints(monomer_type, imol);

    if (r.first) {
        std::map<std::string, coot::map_index_t> name_map;
        for (int i = 0; i < static_cast<int>(atom_selection.size()); i++) {
            std::string atom_name = coot::atom_id_mmdb_expand(atom_selection[i]->name);
            // In gemmi, altloc is '\0' when there is no alternate location
            // In MMDB, altLoc returns "" (empty string) for no alternate location
            // Original condition: atom_alt_conf.empty() || atom_alt_conf == alt_conf
            // When alt_conf is empty, alt_conf[0] is '\0', so atom->altloc == alt_conf[0]
            // matches atoms with no altloc — equivalent to the empty-string check.
            if (atom_selection[i]->altloc == '\0' || atom_selection[i]->altloc == alt_conf[0]) {
                name_map[atom_name] = coot::map_index_t(i);
            }
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
