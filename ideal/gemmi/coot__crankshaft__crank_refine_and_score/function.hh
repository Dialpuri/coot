#pragma once

// Gemmi port of coot::crankshaft::crank_refine_and_score
// Returns gemmi::Models instead of mmdb::Managers

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <vector>

#include "coot/geometry/residue-and-atom-specs.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_previous_residue/gemmi/function.hh"

namespace coot {
namespace crankshaft {

// Gemmi port: takes gemmi::Structure instead of mmdb::Manager*,
// returns vector of gemmi::Model (by value) instead of mmdb::Manager* pointers.
//
// The crankshaft core algorithms (crankshaft, find_maxima, protein_geometry
// refinement) have no gemmi ports yet, so this port delegates to the existing
// _gemmi variants where available and returns empty results when the core
// crankshaft logic cannot be invoked without mmdb types.
std::vector<gemmi::Model>
crank_refine_and_score_gemmi(
  const coot::residue_spec_t &rs,
  unsigned int n_peptides,
  const clipper::Xmap<float> &xmap,
  const gemmi::Structure &st,
  float map_weight,
  int n_samples,
  int n_solutions,
  int n_threads = 0
) {
  std::vector<gemmi::Model> solution_molecules;

  gemmi::Residue *prev_res = coot::util::get_previous_residue_gemmi(rs, st);
  if (!prev_res) {
    // No previous residue — can't run crankshaft
    return solution_molecules;
  }

  // The crankshaft core (find_maxima, protein_geometry refinement,
  // molecule scoring) requires mmdb types and has no gemmi port yet.
  // Until those are ported, return empty solutions.
  // (The original function would call cs.find_maxima(), create models
  // with cs.new_mol_with_moved_atoms(), refine them, score them, and
  // return the top n_solutions.)

  return solution_molecules;
}

} // namespace crankshaft
} // namespace coot