#include "function.hh"
#include <iomanip>
#include <ios>

std::ostream& print_mat33_gemmi(std::ostream& c, gemmi::Mat33 a) {
  c.flags(std::ios::fixed | std::ios::right);
  c.precision(6);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      c << std::setw(10) << a.a[i][j] << " ";
    }
    c << "\n";
  }

  return c;
}