#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "geometry/protein-geometry.hh"
#include "coot-utils/bonded-pairs.hh"
#include "coot/utils/coot-utils.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__find_peptide_link_type_ng/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__make_rama_triples/gemmi/function.hh"

namespace coot {

/// gemmi-port of coot::restraints_container_t::make_rama_plot_restraints_ng
///
/// Original: iterated residues_vec, checked for consecutive standard amino acids
/// with peptide-bond links (TRANS/PTRANS), and called add_rama().
///
/// Port: takes gemmi CRA list (carries chain context) and fixed flags,
/// returns a vector of rama_triple_t_gemmi for the caller to consume.
inline std::vector<rama_triple_t_gemmi> make_rama_plot_restraints_ng_gemmi(
    const std::vector<gemmi::CRA>& residues,
    const std::vector<bool>& fixed_flags,
    const std::map<gemmi::Residue*, std::vector<gemmi::Residue*>>& residue_link_vector_map,
    const std::set<std::pair<gemmi::Residue*, gemmi::Residue*>>& residue_pair_link_set,
    const protein_geometry& geom)
{
    std::vector<rama_triple_t_gemmi> triples;

    int n = static_cast<int>(residues.size());
    if (n < 3) return triples;

    // Mirror the original loop: i-1 (prev), i (this), i+1 (next)
    // Original used index_p/t/n derived from residue->index;
    // with a sorted residues_vec the vector positions ARE the indices.
    for (int i = 1; i < n - 1; ++i) {
        const gemmi::CRA& cra_prev = residues[i - 1];
        const gemmi::CRA& cra_this = residues[i];
        const gemmi::CRA& cra_next = residues[i + 1];

        gemmi::Residue* residue_prev_p = cra_prev.residue;
        gemmi::Residue* residue_this_p = cra_this.residue;
        gemmi::Residue* residue_next_p = cra_next.residue;

        if (!residue_prev_p || !residue_this_p || !residue_next_p)
            continue;

        bool f1 = fixed_flags[i - 1];
        bool f2 = fixed_flags[i];
        bool f3 = fixed_flags[i + 1];

        // If all three are fixed, skip
        if (f1 && f2 && f3)
            continue;

        // Consecutive positions ⇒ consecutive residues in the chain
        // (original checked residue->index differences == 1)
        int index_p = i - 1;
        int index_t = i;
        int index_n = i + 1;

        if ((index_t - index_p) != 1) continue;
        if ((index_n - index_t) != 1) continue;

        // All three must be standard amino acids
        std::string rn_1 = residue_prev_p->name;
        std::string rn_2 = residue_this_p->name;
        std::string rn_3 = residue_next_p->name;

        if (!coot::util::is_standard_amino_acid_name(rn_1)) continue;
        if (!coot::util::is_standard_amino_acid_name(rn_2)) continue;
        if (!coot::util::is_standard_amino_acid_name(rn_3)) continue;

        // Check peptide-bond links via the map
        auto it_map = residue_link_vector_map.find(residue_this_p);
        if (it_map == residue_link_vector_map.end()) continue;

        const std::vector<gemmi::Residue*>& neighbors = it_map->second;

        auto prev_it = std::find(neighbors.begin(), neighbors.end(), residue_prev_p);
        auto next_it = std::find(neighbors.begin(), neighbors.end(), residue_next_p);

        if (prev_it == neighbors.end() || next_it == neighbors.end()) continue;

        // Determine link type — only TRANS and PTRANS add rama restraints
        std::string link_type = coot::restraints_container_t::find_peptide_link_type_ng_gemmi(
            *residue_prev_p, *residue_this_p, geom);

        if (link_type == "TRANS" || link_type == "PTRANS") {
            triples.emplace_back(
                cra_prev.chain,  // chain_1 (from CRA — gemmi Residue has no parent pointer)
                residue_prev_p,
                residue_this_p,
                residue_next_p,
                link_type,
                f1, f2, f3
            );
        }
    }

    return triples;
}

} // namespace coot
