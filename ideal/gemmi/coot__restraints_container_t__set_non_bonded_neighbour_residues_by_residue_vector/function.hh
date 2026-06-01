#pragma once

#include <gemmi/model.hpp>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "geometry/protein-geometry.hh"

namespace coot {

// gemmi-compatible bonded pair (replaces coot::bonded_pair_t which uses mmdb::Residue*)
struct bonded_pair_gemmi {
    gemmi::Residue *res_1;
    gemmi::Residue *res_2;
    bool is_fixed_first;
    bool is_fixed_second;
    std::string link_type;

    bonded_pair_gemmi(gemmi::Residue *r1,
                      gemmi::Residue *r2,
                      bool f1,
                      bool f2,
                      const std::string &lt)
        : res_1(r1), res_2(r2), is_fixed_first(f1), is_fixed_second(f2), link_type(lt) {}
};

// gemmi port of coot::restraints_container_t::set_non_bonded_neighbour_residues_by_residue_vector
//
// Original: member method that sets restraints_container_t::non_bonded_neighbour_residues.
// Ported:   free function that returns the resulting vector (gemmi::Residue* instead of
//           mmdb::Residue*, no mmdb dependency).
inline std::vector<gemmi::Residue *>
set_non_bonded_neighbour_residues_by_residue_vector_gemmi(
    const std::map<gemmi::Residue *, std::set<gemmi::Residue *>> &neighbour_set,
    const std::vector<std::pair<bool, gemmi::Residue *>> &residues_vec,
    const std::vector<bonded_pair_gemmi> &bonded_flanking_pairs,
    const coot::protein_geometry &geom)
{
    (void)geom; // unused (same as original)

    std::vector<gemmi::Residue *> nbr; // non-bonded residues

    for (auto it_map = neighbour_set.begin(); it_map != neighbour_set.end(); ++it_map) {

        const std::set<gemmi::Residue *> &neighbours = it_map->second;
        std::set<gemmi::Residue *>::const_iterator it_set;

        for (it_set = neighbours.begin(); it_set != neighbours.end(); ++it_set) {
            gemmi::Residue *test_res = *it_set;
            if (std::find(nbr.begin(), nbr.end(), test_res) == nbr.end()) {
                // not already there...
                bool found = false;

                for (unsigned int ires = 0; ires < residues_vec.size(); ires++) {
                    if (test_res == residues_vec[ires].second) {
                        if (!residues_vec[ires].first) {
                            found = true;
                            break;
                        }
                    }
                }

                if (!found) {
                    // OK, so this neighbour was not in the passed set of
                    // moving residues (and not already in nbr)... it can
                    // be a flanking residue then...

                    // check that it is not a bonded flanking residue...
                    for (unsigned int iflank = 0; iflank < bonded_flanking_pairs.size(); iflank++) {
                        if (bonded_flanking_pairs[iflank].res_1 == test_res) {
                            found = true;
                            break;
                        }
                        if (bonded_flanking_pairs[iflank].res_2 == test_res) {
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        nbr.push_back(test_res);
                    }
                }
            }
        }
    }
    return nbr;
}

} // namespace coot