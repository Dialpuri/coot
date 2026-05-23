#pragma once

#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/polyheur.hpp>

namespace coot { namespace util {

inline
std::pair<gemmi::Chain*, gemmi::Structure>
copy_chain_gemmi(const gemmi::Chain& chain) {
    gemmi::Structure st;
    st.models.resize(1);
    gemmi::Model& model = st.models[0];
    
    gemmi::Chain new_chain = chain;
    new_chain.name = chain.name;  // copy chain ID
    new_chain.residues = chain.residues;  // copy all residues
    
    model.chains.push_back(new_chain);
    
    gemmi::setup_entities(st);
    
    return std::make_pair(&model.chains.back(), std::move(st));
}

}} // namespace coot::util