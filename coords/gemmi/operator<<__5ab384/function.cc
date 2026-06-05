#include "function.hh"

std::ostream& operator_stream_insert_gemmi(std::ostream& s, SymmMatrix m) {
   s << m.mat[0][0] << " "  << m.mat[0][1] << " "
     << m.mat[0][2] << " "  << m.mat[0][3] << " "  << std::endl
     << m.mat[1][0] << " "  << m.mat[1][1] << " "
     << m.mat[1][2] << " "  << m.mat[1][3] << " "  << std::endl
     << m.mat[2][0] << " "  << m.mat[2][1] << " "
     << m.mat[2][2] << " "  << m.mat[2][3] << " "  << std::endl
     << m.mat[3][0] << " "  << m.mat[3][1] << " "
     << m.mat[3][2] << " "  << m.mat[3][3] << " "  << std::endl;
   return s;
}