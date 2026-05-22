#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Returns the nth residue in the structure (1-indexed), traversing chains in order.
// Returns a CRA with nullptr residue if nth <= 0 or if there are fewer than nth residues.
inline gemmi::CRA get_nth_residue_gemmi(int nth, const gemmi::Structure& st) {
    int count = 0;
    gemmi::CRA result{nullptr, nullptr, nullptr};
    
    if (!st.models.empty()) {
        const gemmi::Model& model = st.models[0];
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                count++;
                if (count == nth) {
                    result = {const_cast<gemmi::Chain*>(&chain),
                              const_cast<gemmi::Residue*>(&residue),
                              nullptr};
                    return result;
                }
            }
        }
    }
    return result;
}

} // namespace util
} // namespace coot