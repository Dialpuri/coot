#pragma once

#include <gemmi/model.hpp>
#include <utility>

namespace coot { namespace util {

inline std::pair<gemmi::Chain*, gemmi::Structure> copy_chain_gemmi(const gemmi::Chain& chain) {
    gemmi::Structure copied_st;
    gemmi::Model copied_model;
    gemmi::Chain copied_chain = chain;  // Deep copy via copy constructor (vectors are copied)
    copied_model.chains.push_back(std::move(copied_chain));
    copied_st.models.push_back(std::move(copied_model));
    gemmi::Chain* copied_chain_ptr = &(copied_st.models.front().chains.front());
    return std::make_pair(copied_chain_ptr, std::move(copied_st));
}

inline void replace_chain_contents_with_atoms_from_chain_gemmi(
    gemmi::Chain* orig_from_chain,
    gemmi::Structure* /*orig_structure*/,
    const gemmi::Chain* modified_chain,
    bool /*do_finishstructedit*/)
{
    if (orig_from_chain && modified_chain) {
        *orig_from_chain = *modified_chain;
    }
}

}} // namespace coot::util