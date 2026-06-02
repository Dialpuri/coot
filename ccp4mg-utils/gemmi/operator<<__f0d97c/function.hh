#pragma once

#include <gemmi/math.hpp>
#include <iosfwd>
#include <ios>
#include <iomanip>

std::ostream& print_mat33_gemmi(std::ostream& c, gemmi::Mat33 a);