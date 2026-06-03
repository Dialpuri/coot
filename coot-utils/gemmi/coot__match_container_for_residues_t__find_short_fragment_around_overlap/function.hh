#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <vector>
#include <string>

namespace coot {

// gemmi port of delete_a_chain_t
struct delete_a_chain_t {
    enum delete_a_chain_type_t {
        NONE = 0,
        DELETE_FIRST_CHAIN = 1,
        DELETE_SECOND_CHAIN = 2
    };

    bool chains_were_mergeable;
    bool short_fragment_is_in_first_selection;
    bool short_fragment_is_upstream_fragment;
    delete_a_chain_type_t delete_type;

    delete_a_chain_t()
        : chains_were_mergeable(false),
          short_fragment_is_in_first_selection(false),
          short_fragment_is_upstream_fragment(false),
          delete_type(NONE) {}

    delete_a_chain_t(bool mergeable, bool in_first, bool upstream)
        : chains_were_mergeable(mergeable),
          short_fragment_is_in_first_selection(in_first),
          short_fragment_is_upstream_fragment(upstream),
          delete_type(NONE) {}
};

// gemmi port of match_container_for_residues_t
class match_container_for_residues_t {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> atom_pairs;

    match_container_for_residues_t() : residue_1(nullptr), residue_2(nullptr) {}

    delete_a_chain_t find_short_fragment_around_overlap_gemmi(
        const std::vector<gemmi::Atom*>& atom_selection_1,
        const std::vector<gemmi::Atom*>& atom_selection_2) const;
};

} // namespace coot