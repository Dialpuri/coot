#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "ideal/parallel-planes.hh"

namespace coot {

inline std::ostream& stream_parallel_planes_t_gemmi(std::ostream& s, coot::parallel_planes_t pp) {
   s << "pp-restr: " << pp.plane_1_atoms.res_spec << " " << pp.plane_2_atoms.res_spec;
   return s;
}

} // namespace coot