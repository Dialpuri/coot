#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

std::vector<gemmi::Atom*> atoms_with_zero_occupancy_gemmi(gemmi::Model& model);

inline std::vector<gemmi::Atom*> atoms_with_zero_occupancy_gemmi(gemmi::Model& model) {
    std::vector<gemmi::Atom*> v;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& at : res.atoms) {
                // In gemmi, TER records are not stored as atoms — they are implicit
                // chain boundaries. So there's no equivalent of isTer() to check.
                if (at.occ < 0.01) {
                    v.push_back(&at);
                }
            }
        }
    }
    return v;
}

} // namespace coot