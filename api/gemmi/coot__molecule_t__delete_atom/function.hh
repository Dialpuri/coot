#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/calculate.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

std::pair<int, size_t> delete_atom_gemmi(gemmi::Structure& st, const atom_spec_t& atom_spec);

}