#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <map>

#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residue_spec_t__select_atoms/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

/* Gemmi port of restraints_container_t::make_h_bond_restraints_from_res_vec_auto.
 *
 * Takes residue info (SeqId, resname, chain_id) instead of mmdb::Residue pointers,
 * finds H-bonds using gemmi::NeighborSearch, and returns results.
 */
struct h_bond_result {
    gemmi::Atom* donor;
    gemmi::Atom* acceptor;
    double dist;
};

std::vector<h_bond_result>
make_h_bond_restraints_from_res_vec_auto_gemmi(
    const std::vector<gemmi::SeqId>& seq_ids,
    const std::vector<std::string>& res_names,
    const std::vector<std::string>& chain_ids,
    gemmi::Structure* st,
    const coot::protein_geometry& geom);

} // namespace coot
