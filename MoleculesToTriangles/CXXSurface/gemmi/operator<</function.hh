#pragma once
#include <ostream>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

namespace coot {

template <typename T>
std::ostream& coot_operator_gemmi(std::ostream& out, const CXXCoord<T>& c) {
    out << "[ ";
    for (int k = 0; k < 4; k++) {
        out << c[k] << ' ';
    }
    out << ']';
    return out;
}

} // namespace coot