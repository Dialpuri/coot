#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <utility>
#include <vector>

#include "geometry/hb-types.hh"
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__h_bonds__mark_donors_and_acceptors/gemmi/function.hh"

// gemmi port of coot::h_bonds::check_hb_status
//
// Original (MMDB):
//   std::pair<bool, int> check_hb_status(int selhnd, mmdb::Manager* mol,
//                                         const protein_geometry& geom, int imol)
//
// The MMDB version used a UDD handle (int) to store per-atom H-bond types.
// The gemmi version returns the map directly (from mark_donors_and_acceptors_gemmi),
// so the second element is the map instead of a handle integer.
//
// sel_atoms  — atoms in the selection (CRA for context)
// geom       — protein geometry with H-bond type tables
// imol       — molecule index (passed through to geometry lookups)
//
// Returns {bool, map}:
//   .first  — true if any atom was marked as donor/acceptor (not HB_UNASSIGNED)
//   .second — the full per-atom H-bond type map from mark_donors_and_acceptors_gemmi

namespace coot {
namespace h_bonds {

inline std::pair<bool, std::map<gemmi::Atom*, int>> check_hb_status_gemmi(
    const std::vector<gemmi::CRA>& sel_atoms,
    const coot::protein_geometry& geom,
    int imol)
{
    bool status = false;

    // Call the gemmi port of mark_donors_and_acceptors with an empty
    // second selection (equivalent to passing -1 as selhnd in the original).
    std::map<gemmi::Atom*, int> hb_map = mark_donors_and_acceptors_gemmi(
        sel_atoms, {}, geom, imol);

    // Check if any atom was assigned a donor/acceptor type
    for (const auto& kv : hb_map) {
        if (kv.second != HB_UNASSIGNED) {
            status = true;
            break;
        }
    }

    return {status, hb_map};
}

} // namespace h_bonds
} // namespace coot