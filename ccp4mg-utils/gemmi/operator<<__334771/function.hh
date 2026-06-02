#pragma once

#include "ccp4mg-utils/quat.h"
#include <iostream>

inline std::ostream& print_quat_gemmi(std::ostream& c, Quat a) {
  const double* v = a.Getdval();
  c << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
  return c;
}