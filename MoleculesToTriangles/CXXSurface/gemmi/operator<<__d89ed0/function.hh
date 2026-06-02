#pragma once

#include <iostream>
#include <sstream>

template <typename T>
class CXXCoord {
public:
    T xyzr[4];

    CXXCoord() {
        for (int i = 0; i < 4; i++) xyzr[i] = 0.;
    }
    CXXCoord(const T *x) {
        for (int i = 0; i < 4; i++) xyzr[i] = x[i];
    }
    CXXCoord(const T x, const T y, const T z) {
        xyzr[0] = x; xyzr[1] = y; xyzr[2] = z; xyzr[3] = 0.;
    }
    CXXCoord(const T x, const T y, const T z, const T r) {
        xyzr[0] = x; xyzr[1] = y; xyzr[2] = z; xyzr[3] = r;
    }
    const T& operator[] (unsigned element) const { return xyzr[element]; }
    T& operator[] (unsigned element) { return xyzr[element]; }
};

template <typename T>
std::ostream& stream_CXXCoord_gemmi(std::ostream &out, const CXXCoord<T>& c)
{
    out << "[ ";
    for (int k = 0; k < 4; k++) {
        out << c[k] << ' ';
    }
    out << ']';
    return out;
}