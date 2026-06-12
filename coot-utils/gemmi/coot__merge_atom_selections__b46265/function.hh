#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/to_pdb.hpp>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

// ── Verified _gemmi port dependencies (minimal non-conflicting set) ────
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__mergeable_atom_selections/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__delete_the_matched_residues_matched_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__renumber_chains_start_at_least_at_1/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__pdbcleanup_serial_residue_numbers/gemmi/function.hh"

namespace coot {

// ── delete_a_chain_t (inline from coot original) ──────────────────────

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

    delete_a_chain_t(const delete_a_chain_t& other)
        : chains_were_mergeable(other.chains_were_mergeable),
          short_fragment_is_in_first_selection(other.short_fragment_is_in_first_selection),
          short_fragment_is_upstream_fragment(other.short_fragment_is_upstream_fragment),
          delete_type(other.delete_type) {}

    delete_a_chain_t& operator=(const delete_a_chain_t& other) {
        if (this != &other) {
            chains_were_mergeable = other.chains_were_mergeable;
            short_fragment_is_in_first_selection = other.short_fragment_is_in_first_selection;
            short_fragment_is_upstream_fragment = other.short_fragment_is_upstream_fragment;
            delete_type = other.delete_type;
        }
        return *this;
    }
};

// ── Helper: find short fragment around overlap ─────────────────────────

static delete_a_chain_t find_short_fragment_around_overlap(
    const match_container_for_residues_t_gemmi& m,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2) {

    delete_a_chain_t result;

    if (m.residue_1 == nullptr || m.residue_2 == nullptr) {
        result.chains_were_mergeable = false;
        return result;
    }

    auto count_residues = [](const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel) -> size_t {
        std::set<const gemmi::Residue*> seen;
        for (auto& p : sel) {
            if (p.second) seen.insert(p.second);
        }
        return seen.size();
    };

    size_t n1 = count_residues(sel1);
    size_t n2 = count_residues(sel2);

    result.chains_were_mergeable = true;
    result.short_fragment_is_in_first_selection = (n1 <= n2);

    int seq1 = m.residue_1 ? m.residue_1->seqid.num.value : 0;
    int seq2 = m.residue_2 ? m.residue_2->seqid.num.value : 0;
    result.short_fragment_is_upstream_fragment = (seq1 < seq2);

    // Determine delete_type
    if (result.short_fragment_is_in_first_selection) {
        result.delete_type = delete_a_chain_t::DELETE_FIRST_CHAIN;
    } else {
        result.delete_type = delete_a_chain_t::DELETE_SECOND_CHAIN;
    }

    return result;
}

// ── Helper: meld chains ────────────────────────────────────────────────

static void meld_chains_gemmi(
    gemmi::Model& model,
    gemmi::Residue* residue_1,
    gemmi::Residue* residue_2,
    bool short_in_first) {

    if (!residue_1 || !residue_2) return;

    gemmi::Chain* chain1 = nullptr;
    gemmi::Chain* chain2 = nullptr;
    for (auto& ch : model.chains) {
        for (auto& res : ch.residues) {
            if (&res == residue_1) chain1 = &ch;
            if (&res == residue_2) chain2 = &ch;
        }
        if (chain1 && chain2) break;
    }
    if (!chain1 || !chain2) return;
    if (chain1 == chain2) return;

    gemmi::Chain* keep_chain   = short_in_first ? chain2 : chain1;
    gemmi::Chain* absorb_chain = short_in_first ? chain1 : chain2;
    gemmi::Residue* overlap_res = short_in_first ? residue_1 : residue_2;

    // Move non-overlap residues from absorb_chain to keep_chain
    std::vector<gemmi::Residue> to_move;
    for (auto it = absorb_chain->residues.begin(); it != absorb_chain->residues.end(); ) {
        if (&(*it) == overlap_res) {
            ++it;
        } else {
            to_move.push_back(std::move(*it));
            it = absorb_chain->residues.erase(it);
        }
    }

    for (auto& res : to_move) {
        keep_chain->residues.push_back(std::move(res));
    }

    // Sort by sequence number
    std::sort(keep_chain->residues.begin(), keep_chain->residues.end(),
              [](const gemmi::Residue& a, const gemmi::Residue& b) {
                  if (a.seqid.num.value != b.seqid.num.value)
                      return a.seqid.num.value < b.seqid.num.value;
                  return a.seqid.icode < b.seqid.icode;
              });

    // Remove empty absorb_chain
    auto cit = std::find_if(model.chains.begin(), model.chains.end(),
        [absorb_chain](const gemmi::Chain& ch) { return &ch == absorb_chain; });
    if (cit != model.chains.end() && cit->residues.empty()) {
        model.chains.erase(cit);
    }
}

// ── Helper: delete residues upstream from matcher ──────────────────────

static void delete_upstream_gemmi(
    gemmi::Model& model,
    gemmi::Residue* matcher_residue,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel_atoms) {

    if (!matcher_residue) return;

    gemmi::Chain* chain_ptr = nullptr;
    for (auto& ch : model.chains) {
        for (auto& res : ch.residues) {
            if (&res == matcher_residue) {
                chain_ptr = &ch;
                break;
            }
        }
        if (chain_ptr) break;
    }
    if (!chain_ptr) return;

    // Build set of residues to delete
    std::set<gemmi::Residue*> to_delete;
    for (auto& p : sel_atoms) {
        if (p.second && p.second != matcher_residue) {
            // Check atom is in this chain
            for (auto& res : chain_ptr->residues) {
                for (auto& atom : res.atoms) {
                    if (&atom == p.first) {
                        to_delete.insert(p.second);
                        break;
                    }
                }
            }
        }
    }

    if (!to_delete.empty()) {
        chain_ptr->residues.erase(
            std::remove_if(chain_ptr->residues.begin(), chain_ptr->residues.end(),
                [&to_delete](const gemmi::Residue& r) {
                    return to_delete.find(&const_cast<gemmi::Residue&>(r)) != to_delete.end();
                }),
            chain_ptr->residues.end());
    }
}

// ── Helper: delete residues downstream from matcher ────────────────────

static void delete_downstream_gemmi(
    gemmi::Model& model,
    gemmi::Residue* matcher_residue,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel_atoms) {

    if (!matcher_residue) return;

    gemmi::Chain* chain_ptr = nullptr;
    for (auto& ch : model.chains) {
        for (auto& res : ch.residues) {
            if (&res == matcher_residue) {
                chain_ptr = &ch;
                break;
            }
        }
        if (chain_ptr) break;
    }
    if (!chain_ptr) return;

    std::set<gemmi::Residue*> to_delete;
    for (auto& p : sel_atoms) {
        if (p.second && p.second != matcher_residue) {
            to_delete.insert(p.second);
        }
    }

    if (!to_delete.empty()) {
        chain_ptr->residues.erase(
            std::remove_if(chain_ptr->residues.begin(), chain_ptr->residues.end(),
                [&to_delete](const gemmi::Residue& r) {
                    return to_delete.find(&const_cast<gemmi::Residue&>(r)) != to_delete.end();
                }),
            chain_ptr->residues.end());
    }
}

// ── Helper: delete short overlapping chain ─────────────────────────────

static void delete_short_overlapping_chain(
    const delete_a_chain_t& merge_values,
    gemmi::Structure& st,
    const std::string& chain_id_i,
    const std::string& chain_id_j) {

    if (st.models.empty()) return;
    gemmi::Model& model = st.models[0];

    std::string chain_to_delete;
    if (merge_values.delete_type == delete_a_chain_t::DELETE_FIRST_CHAIN) {
        chain_to_delete = chain_id_i;
    } else if (merge_values.delete_type == delete_a_chain_t::DELETE_SECOND_CHAIN) {
        chain_to_delete = chain_id_j;
    } else {
        return;
    }

    auto it = std::find_if(model.chains.begin(), model.chains.end(),
        [&chain_to_delete](const gemmi::Chain& ch) { return ch.name == chain_to_delete; });
    if (it != model.chains.end()) {
        model.chains.erase(it);
    }
}

// ── Top-level merge_atom_selections_gemmi ──────────────────────────────

inline void merge_atom_selections_gemmi(gemmi::Structure& st) {
    if (st.models.empty()) return;

    gemmi::Model& model = st.models[0];

    bool continue_merging = true;
    std::set<std::string> merged_chain_ids;

    std::cout << "................ start looping" << std::endl;

    while (continue_merging) {
        std::cout << "................ new loop" << std::endl;
        continue_merging = false;

        // Build selections: one per chain
        std::vector<std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>> selections;
        std::vector<std::string> chain_ids;

        for (gemmi::Chain& chain : model.chains) {
            chain_ids.push_back(chain.name);
            std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>> sel;
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    sel.push_back({&atom, &res});
                }
            }
            selections.push_back(sel);
        }

        bool r = false;

        for (unsigned int i = 0; i < selections.size(); i++) {
            for (unsigned int j = 0; j < selections.size(); j++) {
                if (i != j) {
                    std::string chain_id_chain_id_key = chain_ids[i] + "-" + chain_ids[j];
                    if (merged_chain_ids.find(chain_id_chain_id_key) != merged_chain_ids.end())
                        continue;

                    // Check if mergeable
                    auto m = coot::mergeable_atom_selections_gemmi(selections[i], selections[j]);

                    if (!m.first) continue;

                    // Find short fragment info
                    delete_a_chain_t merge_values = find_short_fragment_around_overlap(
                        m.second, selections[i], selections[j]);

                    if (merge_values.chains_were_mergeable) {
                        r = true;

                        // Delete the matched residues
                        coot::delete_the_matched_residues_matched_residue_gemmi(
                            st, m.second.residue_1, m.second.residue_2,
                            merge_values.short_fragment_is_in_first_selection);

                        // Meld the chains
                        meld_chains_gemmi(st.models[0], m.second.residue_1, m.second.residue_2,
                                        merge_values.short_fragment_is_in_first_selection);

                        // Delete short fragment residues
                        if (merge_values.short_fragment_is_in_first_selection) {
                            if (merge_values.short_fragment_is_upstream_fragment) {
                                delete_upstream_gemmi(st.models[0], m.second.residue_1, selections[i]);
                            } else {
                                delete_downstream_gemmi(st.models[0], m.second.residue_1, selections[i]);
                            }
                        } else {
                            if (merge_values.short_fragment_is_upstream_fragment) {
                                delete_upstream_gemmi(st.models[0], m.second.residue_2, selections[j]);
                            } else {
                                delete_downstream_gemmi(st.models[0], m.second.residue_2, selections[j]);
                            }
                        }

                        // Delete short chain
                        delete_short_overlapping_chain(merge_values, st, chain_ids[i], chain_ids[j]);

                        // Renumber and clean up
                        coot::renumber_chains_start_at_least_at_1_gemmi(st);
                        coot::util::pdbcleanup_serial_residue_numbers_gemmi(st);

                        merged_chain_ids.insert(chain_id_chain_id_key);
                        continue_merging = true;
                        break;
                    }
                }
            }
            if (r) break;
        }
    }
}

} // namespace coot
