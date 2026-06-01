#pragma once

#include <gemmi/model.hpp>
#include <gemmi/seqid.hpp>
#include <utility>
#include "coot/geometry/residue-and-atom-specs.hh"

namespace coot {

std::pair<coot::atom_spec_t, coot::atom_spec_t>
link_atoms_gemmi(const gemmi::Connection& link, gemmi::Model* model_p = nullptr);

} // namespace coot