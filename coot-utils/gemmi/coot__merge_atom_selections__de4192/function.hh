#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <limits>

namespace coot {

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
        : chains_were_mergeable(false), short_fragment_is_in_first_selection(false),
          short_fragment_is_upstream_fragment(false), delete_type(NONE) {}

    delete_a_chain_t(bool mergeable, bool in_first, bool upstream)
        : chains_were_mergeable(mergeable),
          short_fragment_is_in_first_selection(in_first),
          short_fragment_is_upstream_fragment(upstream),
          delete_type(NONE) {}
};

struct match_container_for_residues_atom_pair_t_gemmi {
    gemmi::CRA first;
    gemmi::CRA second;
};

struct match_container_for_residues_t_gemmi {
    gemmi::Residue* residue_1;
    gemmi::Residue* residue_2;
    std::vector<match_container_for_residues_atom_pair_t_gemmi> atom_pairs;

    match_container_for_residues_t_gemmi()
        : residue_1(nullptr), residue_2(nullptr) {}
};

// ─── Helper: compute residue center ───────────────────────────────────

inline gemmi::Vec3 _residue_center(const gemmi::Residue& res) {
    if (res.atoms.empty()) return gemmi::Vec3{0, 0, 0};
    gemmi::Vec3 sum{0, 0, 0};
    for (const auto& atom : res.atoms) {
        sum += atom.pos;
    }
    return sum / static_cast<double>(res.atoms.size());
}

// ─── mergeable_atom_selections_gemmi ──────────────────────────────────

inline std::pair<bool, match_container_for_residues_t_gemmi>
mergeable_atom_selections_gemmi(
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2) {

    match_container_for_residues_t_gemmi mc;

    if (sel1.empty() || sel2.empty())
        return {false, mc};

    std::map<std::string, gemmi::Atom*> atoms1;
    std::map<std::string, gemmi::Atom*> atoms2;

    for (const auto& ap : sel1)
        if (ap.first) atoms1[ap.first->name] = ap.first;
    for (const auto& ap : sel2)
        if (ap.first) atoms2[ap.first->name] = ap.first;

    std::vector<match_container_for_residues_atom_pair_t_gemmi> matches;
    for (auto& [name, a1] : atoms1) {
        auto it = atoms2.find(name);
        if (it != atoms2.end()) {
            double dist_sq = (a1->pos - it->second->pos).length_sq();
            if (dist_sq < 0.02 * 0.02) {
                match_container_for_residues_atom_pair_t_gemmi apair;
                for (const auto& sp : sel1) {
                    if (sp.first == a1) {
                        apair.first = gemmi::CRA{nullptr, sp.second, a1};
                        break;
                    }
                }
                for (const auto& sp : sel2) {
                    if (sp.first == it->second) {
                        apair.second = gemmi::CRA{nullptr, sp.second, it->second};
                        break;
                    }
                }
                matches.push_back(apair);
            }
        }
    }

    if (matches.empty())
        return {false, mc};

    std::map<std::pair<gemmi::Residue*, gemmi::Residue*>, int> residue_match_count;
    for (const auto& m : matches) {
        residue_match_count[{m.first.residue, m.second.residue}]++;
    }

    gemmi::Residue* r1 = nullptr;
    gemmi::Residue* r2 = nullptr;
    int best_count = 0;
    for (const auto& [key, count] : residue_match_count) {
        if (count > best_count) {
            best_count = count;
            r1 = key.first;
            r2 = key.second;
        }
    }

    if (!r1 || !r2)
        return {false, mc};

    for (const auto& m : matches) {
        if (m.first.residue != r1 || m.second.residue != r2)
            return {false, mc};
    }

    double dist = _residue_center(*r1).dist(_residue_center(*r2));
    if (dist > 0.8)
        return {false, mc};

    mc.residue_1 = r1;
    mc.residue_2 = r2;
    mc.atom_pairs = matches;
    return {true, mc};
}

// ─── find_short_fragment_around_overlap_gemmi ─────────────────────────

inline coot::delete_a_chain_t
find_short_fragment_around_overlap_gemmi(
    match_container_for_residues_t_gemmi m,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2) {

    coot::delete_a_chain_t result;

    std::set<gemmi::Residue*> residues1;
    std::set<gemmi::Residue*> residues2;
    for (const auto& ap : sel1) if (ap.second) residues1.insert(ap.second);
    for (const auto& ap : sel2) if (ap.second) residues2.insert(ap.second);

    gemmi::Residue* overlap_res1 = m.residue_1;
    gemmi::Residue* overlap_res2 = m.residue_2;

    int sel1_upstream = 0, sel1_downstream = 0;
    for (gemmi::Residue* r : residues1) {
        if (r == overlap_res1) continue;
        if (r->seqid.num.value < overlap_res1->seqid.num.value)
            sel1_upstream++;
        else
            sel1_downstream++;
    }

    int sel2_upstream = 0, sel2_downstream = 0;
    for (gemmi::Residue* r : residues2) {
        if (r == overlap_res2) continue;
        if (r->seqid.num.value < overlap_res2->seqid.num.value)
            sel2_upstream++;
        else
            sel2_downstream++;
    }

    int sel1_total = static_cast<int>(residues1.size());
    int sel2_total = static_cast<int>(residues2.size());

    if (sel2_total < sel1_total) {
        result.short_fragment_is_in_first_selection = false;
        result.short_fragment_is_upstream_fragment = (sel2_upstream < sel2_downstream);
        result.chains_were_mergeable = true;
    } else if (sel1_total < sel2_total) {
        result.short_fragment_is_in_first_selection = true;
        result.short_fragment_is_upstream_fragment = (sel1_upstream < sel1_downstream);
        result.chains_were_mergeable = true;
    } else {
        if (sel1_total <= 1 && sel2_total <= 1) {
            result.chains_were_mergeable = true;
            result.short_fragment_is_in_first_selection = false;
            result.short_fragment_is_upstream_fragment = true;
        } else {
            result.chains_were_mergeable = false;
        }
    }

    return result;
}

// ─── delete_upstream_gemmi ────────────────────────────────────────────

inline void
delete_upstream_gemmi(
    match_container_for_residues_t_gemmi,
    bool,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>&,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>&) {
    // Selection-level marker; actual deletion happens in meld_gemmi
}

// ─── delete_downstream_gemmi ──────────────────────────────────────────

inline void
delete_downstream_gemmi(
    match_container_for_residues_t_gemmi,
    bool,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>&,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>&) {
    // Selection-level marker; actual deletion happens in meld_gemmi
}

// ─── meld_gemmi ───────────────────────────────────────────────────────

inline void
meld_gemmi(
    gemmi::Structure& st,
    match_container_for_residues_t_gemmi m,
    std::pair<bool, bool> rr) {

    if (st.models.empty()) return;
    gemmi::Model& model = st.models[0];
    if (!m.residue_1 || !m.residue_2) return;

    gemmi::Chain* chain1 = nullptr;
    gemmi::Chain* chain2 = nullptr;
    for (auto& ch : model.chains) {
        for (auto& res : ch.residues) {
            if (&res == m.residue_1) chain1 = &ch;
            if (&res == m.residue_2) chain2 = &ch;
        }
        if (chain1 && chain2) break;
    }
    if (!chain1 || !chain2) return;

    bool short_in_first = rr.first;
    gemmi::Chain* keep_chain   = short_in_first ? chain2 : chain1;
    gemmi::Chain* absorb_chain = short_in_first ? chain1 : chain2;
    gemmi::Residue* overlap_res = short_in_first ? m.residue_1 : m.residue_2;

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

    std::sort(keep_chain->residues.begin(), keep_chain->residues.end(),
              [](const gemmi::Residue& a, const gemmi::Residue& b) {
                  if (a.seqid.num.value != b.seqid.num.value)
                      return a.seqid.num.value < b.seqid.num.value;
                  return a.seqid.icode < b.seqid.icode;
              });
}

// ─── delete_the_matched_residues_matched_residue_gemmi ────────────────

inline void
delete_the_matched_residues_matched_residue_gemmi(
    gemmi::Structure& st,
    gemmi::Residue* residue_1,
    gemmi::Residue* residue_2,
    bool short_fragment_is_in_first_selection) {

    if (st.models.empty()) return;
    gemmi::Model& model = st.models[0];

    gemmi::Residue* keep_res   = short_fragment_is_in_first_selection ? residue_1 : residue_2;
    gemmi::Residue* delete_res = short_fragment_is_in_first_selection ? residue_2 : residue_1;
    if (!delete_res) return;

    for (auto& chain : model.chains) {
        auto it = std::find_if(chain.residues.begin(), chain.residues.end(),
                               [delete_res](const gemmi::Residue& r) {
                                   return &r == delete_res;
                               });
        if (it != chain.residues.end() && &*it != keep_res) {
            chain.residues.erase(it);
            return;
        }
    }
}

// ─── renumber_chains_start_at_least_at_1_gemmi ────────────────────────

inline void
renumber_chains_start_at_least_at_1_gemmi(gemmi::Structure& st) {
    if (st.models.empty()) return;
    gemmi::Model& model = st.models[0];

    for (auto& chain : model.chains) {
        int min_seq = std::numeric_limits<int>::max();
        bool found = false;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value < min_seq) {
                min_seq = res.seqid.num.value;
                found = true;
            }
        }
        if (found && min_seq < 1) {
            int offset = 1 - min_seq;
            for (auto& res : chain.residues) {
                res.seqid.num.value += offset;
            }
        }
    }
}

// ─── merge_atom_selections_gemmi (main function) ──────────────────────

inline coot::delete_a_chain_t merge_atom_selections_gemmi(
    gemmi::Structure& st,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2) {

   coot::delete_a_chain_t merge_values(false, false, false);
   bool done_merge = false;
   std::pair<bool, match_container_for_residues_t_gemmi> m =
       coot::mergeable_atom_selections_gemmi(sel1, sel2);

   if (true)
      std::cout << "DEBUG:: in merge_atom_selections(): "
                << " found mergeable selections status: " << m.first << std::endl;

   if (m.first) {
      coot::delete_a_chain_t r =
          find_short_fragment_around_overlap_gemmi(m.second, sel1, sel2);
      merge_values = r;

      if (true)
         std::cout << "debug:: m.second.find_short_fragment_around_overlap() "
                   << std::endl;

      if (r.chains_were_mergeable) {
         if (r.short_fragment_is_in_first_selection) {
            if (r.short_fragment_is_upstream_fragment) {
               std::cout << "merge_atom_selections(): --- Block A ---" << std::endl;
               delete_upstream_gemmi(m.second, true, sel1, sel2);
               delete_downstream_gemmi(m.second, false, sel1, sel2);
               std::pair<bool,bool> rr(r.short_fragment_is_in_first_selection,
                                       r.short_fragment_is_upstream_fragment);
               meld_gemmi(st, m.second, rr);
            } else {
               std::cout << "merge_atom_selections(): --- Block B ---" << std::endl;
               delete_downstream_gemmi(m.second, true, sel1, sel2);
               delete_upstream_gemmi(m.second, false, sel1, sel2);
               std::pair<bool,bool> rr(r.short_fragment_is_in_first_selection,
                                       r.short_fragment_is_upstream_fragment);
               meld_gemmi(st, m.second, rr);
            }
         } else {
            if (r.short_fragment_is_upstream_fragment) {
               std::cout << "merge_atom_selections(): --- Block C ---" << std::endl;
               delete_upstream_gemmi(m.second, false, sel1, sel2);
               delete_downstream_gemmi(m.second, true, sel1, sel2);
               std::pair<bool,bool> rr(r.short_fragment_is_in_first_selection,
                                       r.short_fragment_is_upstream_fragment);
               meld_gemmi(st, m.second, rr);
            } else {
               std::cout << "merge_atom_selections(): --- Block D ---" << std::endl;
               delete_downstream_gemmi(m.second, false, sel1, sel2);
               delete_upstream_gemmi(m.second, true, sel1, sel2);
               std::pair<bool,bool> rr(r.short_fragment_is_in_first_selection,
                                       r.short_fragment_is_upstream_fragment);
               meld_gemmi(st, m.second, rr);
            }
         }
         done_merge = true;

         std::cout << "debug:: calling delete_the_matched_residues_matched_residue()" << std::endl;
         delete_the_matched_residues_matched_residue_gemmi(
             st, m.second.residue_1, m.second.residue_2,
             r.short_fragment_is_in_first_selection);

         renumber_chains_start_at_least_at_1_gemmi(st);
      }
   }

   return merge_values;
}

} // namespace coot