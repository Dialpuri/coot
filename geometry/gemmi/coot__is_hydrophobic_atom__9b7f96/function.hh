#pragma once

#include <gemmi/model.hpp>
#include "geometry/hydrophobic.hh"

namespace coot {

inline bool is_hydrophobic_atom_gemmi(gemmi::CRA cra) {
    std::string atom_name   = cra.atom ? cra.atom->name   : "";
    std::string res_name    = cra.residue ? cra.residue->name : "";
    return coot::is_hydrophobic_atom(res_name, atom_name);
}

} // namespace coot
