#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Copy a gemmi::Structure by value - returns a deep copy
// This is the gemmi equivalent of the MMDB coot::util::copy_molecule
inline gemmi::Structure copy_molecule_gemmi(const gemmi::Structure& st) {
    return st;
}

}
}