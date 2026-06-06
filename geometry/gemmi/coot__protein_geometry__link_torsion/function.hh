#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "gemmi/cifdoc.hpp"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__link_add_torsion/gemmi/function.hh"

namespace coot::protein_geometry {

void link_torsion_gemmi(const gemmi::cif::Loop &loop);

} // namespace coot::protein_geometry