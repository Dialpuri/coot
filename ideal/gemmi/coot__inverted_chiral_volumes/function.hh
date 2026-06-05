#pragma once
#include <utility>
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <coot/geometry/residue-and-atom-specs.hh>
#include <coot/geometry/protein-geometry.hh>
#include <coot/ideal/simple-restraint.hh>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_inverted_chiral_atom_p/gemmi/function.hh"

namespace coot {

std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t>>
inverted_chiral_volumes_gemmi(int imol,
                              const gemmi::Structure& st,
                              coot::protein_geometry* geom_p,
                              int cif_dictionary_read_number);

} // namespace coot