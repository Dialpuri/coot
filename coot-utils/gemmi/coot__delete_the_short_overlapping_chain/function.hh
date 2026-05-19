#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct delete_a_chain_t {
    enum delete_a_chain_type_t { NONE, DELETE_FIRST_CHAIN, DELETE_SECOND_CHAIN };
    bool chains_were_mergeable;
    bool short_fragment_is_in_first_selection;
    bool short_fragment_is_upstream_fragment;
    delete_a_chain_type_t delete_type;

    delete_a_chain_t(bool a, bool b, bool c)
        : chains_were_mergeable(a), short_fragment_is_in_first_selection(b),
          short_fragment_is_upstream_fragment(c), delete_type(NONE) {}
};

inline void delete_the_short_overlapping_chain_gemmi(delete_a_chain_t dac,
                                                     gemmi::Structure& st,
                                                     const std::string& chain_id_i_chain,
                                                     const std::string& chain_id_j_chain) {
    if (st.models.empty()) return;
    
    gemmi::Model& model = st.models[0];
    std::vector<gemmi::Chain*> chains_to_delete;
    
    for (gemmi::Chain& chain : model.chains) {
        if (dac.delete_type == delete_a_chain_t::DELETE_FIRST_CHAIN) {
            if (chain.name == chain_id_i_chain) {
                chains_to_delete.push_back(&chain);
            }
        }
        if (dac.delete_type == delete_a_chain_t::DELETE_SECOND_CHAIN) {
            if (chain.name == chain_id_j_chain) {
                chains_to_delete.push_back(&chain);
            }
        }
    }
    
    // Delete chains by removing them from the model
    for (gemmi::Chain* chain : chains_to_delete) {
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
            [chain](const gemmi::Chain& c) { return &c == chain; });
        if (it != model.chains.end()) {
            model.chains.erase(it);
        }
    }
}

} // namespace coot