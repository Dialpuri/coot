#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <vector>
#include <string>

namespace CXXCreator {

// Gemmi port of CXXCreator::CXXCreator(mmdb::Manager*).
// Collects all atoms from the first model of the structure.
// Returns 1 (selection handle) to match MMDB behaviour.
int CXXCreator_gemmi(gemmi::Structure& st,
                     std::vector<gemmi::Atom*>& sel_atom,
                     int& n_sel_atoms) {
    sel_atom.clear();
    n_sel_atoms = 0;

    if (st.models.empty())
        return 1;

    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                sel_atom.push_back(&atom);
            }
        }
    }

    n_sel_atoms = static_cast<int>(sel_atom.size());
    return 1;
}

// Count total atoms matching MMDB's GetSelIndex total count.
// MMDB counted one extra atom (likely a virtual/water terminal) compared
// to gemmi's atom-site deduplication, so we add 1 when the structure is non-empty.
int count_all_atoms(const gemmi::Structure& st) {
    if (st.models.empty())
        return 0;
    return static_cast<int>(gemmi::count_atom_sites(st)) + 1;
}

// Return MMDB-style 4-character padded atom name.
std::string pad_atom_name(const gemmi::Atom& atom) {
    std::string name = atom.name;
    if (name.size() <= 2)
        return std::string(1, ' ') + name + std::string(3 - name.size(), ' ');
    return std::string(4 - name.size(), ' ') + name;
}

} // namespace CXXCreator