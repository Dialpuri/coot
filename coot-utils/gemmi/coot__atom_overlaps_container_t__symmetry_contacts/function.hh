#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <algorithm>

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residue_spec_t__select_atoms/gemmi/function.hh"

namespace coot {

struct atom_overlap_t {
  gemmi::Atom* at1_;
  gemmi::Atom* at2_;

  atom_overlap_t(gemmi::Atom* a1, gemmi::Atom* a2)
    : at1_(a1), at2_(a2) {}
};

class atom_overlaps_container_t {
public:
  gemmi::Residue* res_central_;
  std::vector<gemmi::Residue*> neighbour_residues_;
  gemmi::Structure* mol_;
  coot::protein_geometry const* protein_geom_;
  std::string chain_id_;

  atom_overlaps_container_t(
    gemmi::Residue* res,
    std::vector<gemmi::Residue*> neighbours,
    gemmi::Structure* mol,
    coot::protein_geometry const* geom,
    const std::string& chain_id = ""
  )
    : res_central_(res),
      neighbour_residues_(neighbours),
      mol_(mol),
      protein_geom_(geom),
      chain_id_(chain_id)
  {}

  std::vector<atom_overlap_t> symmetry_contacts_gemmi(float d);
};

} // namespace coot
