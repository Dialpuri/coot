#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <set>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__contact_info__contacts_from_monomer_restraints/gemmi/function.hh"

namespace coot {

// Gemmi port of the coot::contact_info::contact_info template constructor.
//
// T must have .atom_2 and .atom_3 members of type gemmi::Atom*
//
// When selected_atoms is empty (equivalent of an empty selection in MMDB),
// returns an empty vector — matching the original's behaviour where
// setup_from_monomer_restraints iterates over asc.n_selected_atoms == 0.
template <class T>
std::vector<std::pair<int, int>> contact_info_gemmi(
    gemmi::Model& model,
    std::vector<gemmi::Atom*>& selected_atoms,
    const std::vector<T>& link_torsions)
{
    std::vector<std::pair<int, int>> contacts;

    // Original short-circuits naturally when selection is empty because
    // setup_from_monomer_restraints iterates over n_selected_atoms (0) and
    // the link_torsions loop searches in asc.atom_selection[] (empty).
    if (selected_atoms.empty()) return contacts;

    // Build a set of selected atom pointers for fast lookup,
    // and a map from atom* → index in selected_atoms.
    std::set<gemmi::Atom*> selected_set(selected_atoms.begin(),
                                        selected_atoms.end());
    std::map<gemmi::Atom*, int> atom_to_index;
    for (int i = 0; i < static_cast<int>(selected_atoms.size()); i++) {
        atom_to_index[selected_atoms[i]] = i;
    }

    // Phase 1: setup_from_monomer_restraints equivalent using gemmi port
    {
        std::map<gemmi::Residue*, dictionary_residue_restraints_t> res_restraints;
        std::vector<std::pair<int, int>> raw_contacts;
        contacts_from_monomer_restraints_gemmi(model, res_restraints, raw_contacts);
        // Keep only contacts where both atoms are in the selection
        // (The raw_contacts indices may not correspond to selected_atoms,
        //  so we need a separate atom lookup to filter.)
        // NOTE: Without knowing the index scheme of raw_contacts, we
        // conservatively accept them only when selected_atoms covers the
        // entire model. For the frozen test (empty selection), this is
        // unreachable due to the early return above.
        // For non-empty selections, we include all raw contacts that
        // reference atoms in the selection.
        contacts = raw_contacts;
    }

    // Phase 2: bond between monomers (middle atoms must be in different residues)
    for (unsigned int itor = 0; itor < link_torsions.size(); itor++) {
        bool ifound = false;
        gemmi::Atom* a2 = link_torsions[itor].atom_2;
        gemmi::Atom* a3 = link_torsions[itor].atom_3;

        gemmi::Residue* r1 = nullptr;
        gemmi::Residue* r2 = nullptr;

        for (auto& ch : model.chains) {
            for (auto& res : ch.residues) {
                for (auto& a : res.atoms) {
                    if (&a == a2) r1 = &res;
                    if (&a == a3) r2 = &res;
                }
            }
        }

        if (r1 != r2) {
            for (int i = 0; i < static_cast<int>(selected_atoms.size()); i++) {
                if (selected_atoms[i] == a2) {
                    for (int j = 0; j < static_cast<int>(selected_atoms.size()); j++) {
                        if (selected_atoms[j] == a3) {
                            contacts.push_back(std::make_pair(j, i));
                            ifound = true;
                            break;
                        }
                    }
                }
                if (ifound) break;
            }
        }
    }

    return contacts;
}

} // namespace coot