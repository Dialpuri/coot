#pragma once

#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include "clipper/core/cell.h"
#include <cmath>

namespace coot {
namespace util {

inline bool set_mol_cell_gemmi(gemmi::Structure& st, clipper::Cell cell_local) {
  // Set the cell: gemmi::UnitCell takes (a,b,c,alpha,beta,gamma) in its constructor.
  // gemmi uses degrees for angles; clipper::Cell stores angles in radians.
  double a = cell_local.a();
  double b = cell_local.b();
  double c = cell_local.c();
  double alpha = clipper::Util::rad2d(cell_local.alpha());
  double beta  = clipper::Util::rad2d(cell_local.beta());
  double gamma = clipper::Util::rad2d(cell_local.gamma());

  st.cell = gemmi::UnitCell(a, b, c, alpha, beta, gamma);

  // Verify the cell was set correctly
  bool status = false;
  if (std::fabs(st.cell.a - cell_local.a()) < 0.1)
    if (std::fabs(st.cell.b - cell_local.b()) < 0.1)
      if (std::fabs(st.cell.c - cell_local.c()) < 0.1)
        if (std::fabs(st.cell.alpha - clipper::Util::rad2d(cell_local.alpha())) < 0.1)
          if (std::fabs(st.cell.beta  - clipper::Util::rad2d(cell_local.beta())) < 0.1)
            if (std::fabs(st.cell.gamma - clipper::Util::rad2d(cell_local.gamma())) < 0.1)
              status = true;

  return status;
}

} // namespace util
} // namespace coot