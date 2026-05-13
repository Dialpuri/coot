#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot { namespace util {

// Reset residue indices so that each residue's index matches its 0-based position within its chain.
// Uses gemmi::Residue::group_idx as the storage location for the index.
inline void pdbcleanup_serial_residue_numbers_gemmi(gemmi::Structure& st) {
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            int nres = static_cast<int>(chain.residues.size());
            for (int ires = 0; ires < nres; ++ires) {
                chain.residues[ires].group_idx = static_cast<short>(ires);
            }
        }
    }
}

}} // namespace coot::util