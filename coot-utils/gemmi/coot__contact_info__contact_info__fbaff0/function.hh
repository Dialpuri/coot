#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <string>
#include <map>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__contact_info__setup_from_monomer_restraints/gemmi/function.hh"

namespace coot {

struct contacts_pair {
    int i;
    int j;
    contacts_pair(int ii, int jj) : i(ii), j(jj) {}
};

struct contact_info {
    std::vector<contacts_pair> contacts;

    contact_info(
        gemmi::Model& model,
        const std::vector<gemmi::Atom*>& atom_selection,
        int imol,
        coot::protein_geometry *geom_p,
        const std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>& link_bond_atoms
    ) {
        // Call the existing gemmi port for monomer restraints setup
        std::vector<std::pair<int, int>> int_contacts;
        setup_from_monomer_restraints_gemmi(model, atom_selection, imol, geom_p, int_contacts);

        // Convert int_contacts to contacts_pair
        for (const auto& p : int_contacts) {
            contacts.push_back(contacts_pair(p.second, p.first));
        }

        // Now the link bond restraints
        for (unsigned int ilb = 0; ilb < link_bond_atoms.size(); ilb++) {
            bool ifound = false;
            for (int i = 0; i < static_cast<int>(atom_selection.size()); i++) {
                if (atom_selection[i] == link_bond_atoms[ilb].first) {
                    for (int j = 0; j < static_cast<int>(atom_selection.size()); j++) {
                        if (atom_selection[j] == link_bond_atoms[ilb].second) {
                            contacts_pair p(j, i);
                            contacts.push_back(p);
                            ifound = true;
                            break;
                        }
                    }
                }
                if (ifound) break;
            }
        }
    }

    int n_contacts() const {
        return static_cast<int>(contacts.size());
    }
};

} // namespace coot