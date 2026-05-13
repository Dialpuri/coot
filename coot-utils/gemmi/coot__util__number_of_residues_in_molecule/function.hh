#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

inline int number_of_residues_in_molecule_gemmi(const gemmi::Structure& st) {
    int number_of_residues = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            number_of_residues += static_cast<int>(chain.residues.size());
        }
    }
    return number_of_residues;
}

}} // namespace coot::util