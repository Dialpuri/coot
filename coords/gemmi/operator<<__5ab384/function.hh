#pragma once

#include <ostream>

/// 4x4 symmetry matrix (mirrors the original MMDB SymmMatrix)
struct SymmMatrix {
  double mat[4][4];
};

/// Gemmi port of operator<<(std::ostream&, SymmMatrix)
/// Prints each row on its own line, values separated by spaces.
std::ostream& operator_stream_insert_gemmi(std::ostream& s, SymmMatrix m);