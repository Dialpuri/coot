#pragma once
#include <gemmi/model.hpp>
#include <vector>

namespace coot { namespace util {

// Port of coot::util::get_hetgroups from MMDB to gemmi
// Returns residues that have het_flag == 'H' (HETATM)
// If include_waters is false, excludes water molecules (HOH)
inline std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>>
get_hetgroups_gemmi(gemmi::Structure& st, bool include_waters) {
    std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>> het_residues;
    
    if (st.models.empty())
        return het_residues;
    
    gemmi::Model& model = st.models[0];
    
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& residue : chain.residues) {
            if (include_waters || !residue.is_water()) {
                if (residue.het_flag == 'H') {
                    het_residues.push_back({&chain, &residue});
                }
            }
        }
    }
    
    return het_residues;
}

}} // namespace coot::util