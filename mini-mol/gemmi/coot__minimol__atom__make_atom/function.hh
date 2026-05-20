#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <coot/mini-mol/mini-mol.hh>

namespace coot {
namespace minimol {

inline gemmi::Atom make_atom_gemmi(const atom &a) {
    gemmi::Atom at;
    at.name = a.name;
    at.element = gemmi::Element(a.element.c_str());
    at.pos.x = a.pos.x();
    at.pos.y = a.pos.y();
    at.pos.z = a.pos.z();
    at.occ = a.occupancy;
    at.b_iso = a.temperature_factor;
    return at;
}

} // namespace minimol
} // namespace coot