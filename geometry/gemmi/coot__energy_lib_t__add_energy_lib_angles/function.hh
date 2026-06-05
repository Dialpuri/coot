#pragma once

#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

namespace coot {

struct energy_lib_angle {
  std::string atom_type_1;
  std::string atom_type_2;
  std::string atom_type_3;
  double spring_constant;
  double angle;
  double angle_esd;

  energy_lib_angle(
    std::string at1,
    std::string at2,
    std::string at3,
    double spring_const,
    double ang,
    double ang_esd
  ) : atom_type_1(at1),
      atom_type_2(at2),
      atom_type_3(at3),
      spring_constant(spring_const),
      angle(ang),
      angle_esd(ang_esd) {}
};

class energy_lib_t {
public:
  std::vector<energy_lib_angle> angles;

  void add_energy_lib_angles_gemmi(const gemmi::cif::Loop& loop);
};

} // namespace coot