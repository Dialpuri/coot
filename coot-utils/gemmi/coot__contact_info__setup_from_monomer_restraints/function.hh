#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "geometry/protein-geometry.hh"

namespace coot {

inline void setup_from_monomer_restraints_gemmi(
    gemmi::Model& model,
    const std::vector<gemmi::Atom*>& atom_selection,
    int imol,
    coot::protein_geometry *geom_p,
    std::vector<std::pair<int, int>>& contacts)
{
    // Build atom→residue mapping (gemmi::Atom has no parent pointer)
    std::map<gemmi::Atom*, gemmi::Residue*> atom_to_residue;
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                atom_to_residue[&atom] = &res;
            }
        }
    }

    // Collect unique residues from the selection
    std::vector<gemmi::Residue*> residues;
    for (int i = 0; i < static_cast<int>(atom_selection.size()); i++) {
        gemmi::Residue* r = atom_to_residue[atom_selection[i]];
        if (std::find(residues.begin(), residues.end(), r) == residues.end()) {
            residues.push_back(r);
        }
    }

    // Get restraints for each unique residue
    std::map<gemmi::Residue*, coot::dictionary_residue_restraints_t> res_restraints;
    for (unsigned int ires = 0; ires < residues.size(); ires++) {
        std::string rn = residues[ires]->name;
        std::pair<bool, coot::dictionary_residue_restraints_t> rest =
            geom_p->get_monomer_restraints(rn, imol);
        if (!rest.first) {
            std::string m = "Restraints not found for type ";
            m += rn;
            throw std::runtime_error(m);
        }
        res_restraints[residues[ires]] = rest.second;
    }

    // Now build contacts: iterate only the selected atoms
    contacts.clear();
    for (size_t iat = 0; iat < atom_selection.size(); iat++) {
        gemmi::Residue* res1 = atom_to_residue[atom_selection[iat]];
        std::string atom_name_1 = coot::atom_id_mmdb_expand(atom_selection[iat]->name);

        for (size_t jat = 0; jat < atom_selection.size(); jat++) {
            gemmi::Residue* res2 = atom_to_residue[atom_selection[jat]];

            // if they are in the same residue...
            if (res1 == res2) {
                std::string atom_name_2 = coot::atom_id_mmdb_expand(atom_selection[jat]->name);

                auto it = res_restraints.find(res1);
                if (it != res_restraints.end()) {
                    const std::vector<dict_bond_restraint_t>& bond_restraints =
                        it->second.bond_restraint;

                    for (size_t ibond = 0; ibond < bond_restraints.size(); ibond++) {
                        if (bond_restraints[ibond].atom_id_1_4c() == atom_name_1) {
                            if (bond_restraints[ibond].atom_id_2_4c() == atom_name_2) {
                                contacts.push_back(std::make_pair(static_cast<int>(iat), static_cast<int>(jat)));
                                break;
                            }
                        }
                        // and the reverse indexing of that...
                        if (bond_restraints[ibond].atom_id_1_4c() == atom_name_2) {
                            if (bond_restraints[ibond].atom_id_2_4c() == atom_name_1) {
                                contacts.push_back(std::make_pair(static_cast<int>(jat), static_cast<int>(iat)));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace coot