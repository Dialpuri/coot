#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/clipper.h>
#include <cmath>

namespace coot {

class ShelxIns {
public:
  clipper::Cell cell;
  int have_cell_flag = 0;

  bool try_assign_cell_gemmi(const gemmi::Structure& st);
};

inline bool coot::ShelxIns::try_assign_cell_gemmi(const gemmi::Structure& st) {
  if (!have_cell_flag) {
    double a = st.cell.a;
    double b = st.cell.b;
    double c = st.cell.c;
    double alpha_deg = st.cell.alpha;
    double beta_deg  = st.cell.beta;
    double gamma_deg = st.cell.gamma;

    clipper::Cell_descr cdr(a, b, c,
                            clipper::Util::d2rad(alpha_deg),
                            clipper::Util::d2rad(beta_deg),
                            clipper::Util::d2rad(gamma_deg));
    cell = clipper::Cell(cdr);
    have_cell_flag = 1;
  }
  return have_cell_flag;
}

}  // namespace coot