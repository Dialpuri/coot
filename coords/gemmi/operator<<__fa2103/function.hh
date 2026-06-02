#pragma once

#include <ostream>

struct Cell_Translation {
    int us, vs, ws;
    Cell_Translation() : us(0), vs(0), ws(0) {}
    Cell_Translation(int a, int b, int c) : us(a), vs(b), ws(c) {}
};

inline std::ostream& stream_cell_translation_gemmi(std::ostream &s, const Cell_Translation &ct) {
    s << "Cell Trans: (" << ct.us << " " << ct.vs << " " << ct.ws << ")";
    return s;
}