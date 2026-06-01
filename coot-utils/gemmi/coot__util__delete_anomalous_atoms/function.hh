#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace util {

/// Remove CB atoms from GLY residues (anomalous: glycine has no beta carbon).
/// Mirrors coot::util::delete_anomalous_atoms but uses gemmi::Structure.
void delete_anomalous_atoms_gemmi(gemmi::Structure& st) {
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.name != "GLY") continue;

                // Collect indices of atoms named " CB " (PDB-space-padded)
                std::vector<int> delete_indices;
                for (int i = 0; i < static_cast<int>(res.atoms.size()); i++) {
                    if (res.atoms[i].name == " CB ") {
                        delete_indices.push_back(i);
                    }
                }

                // Erase in reverse order to preserve indices
                for (auto it = delete_indices.rbegin(); it != delete_indices.rend(); ++it) {
                    res.atoms.erase(res.atoms.begin() + *it);
                }
            }
        }
    }
}

} // namespace util
} // namespace coot