#pragma once

#include "mini-mol/mini-mol.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__residue__residue__7b6771/gemmi/function.hh"

namespace coot {
namespace minimol {

inline molecule molecule_gemmi(const fragment &frag) {
    molecule mol;
    mol.fragments.push_back(frag);
    return mol;
}

} // namespace minimol
} // namespace coot
