#pragma once
#include <iomanip>
#include <ios>
#include <iostream>

struct Cartesian_gemmi {
    double x_, y_, z_, a_;

    Cartesian_gemmi(double x, double y, double z, double a = 1.0)
        : x_(x), y_(y), z_(z), a_(a) {}

    double get_x() const { return x_; }
    double get_y() const { return y_; }
    double get_z() const { return z_; }
    double get_a() const { return a_; }
};

inline std::ostream& operator_stream_gemmi(std::ostream& c, const Cartesian_gemmi& a) {

    c.flags(std::ios::fixed | std::ios::right);
    c.precision(6);

    c << std::setw(10) << a.get_x() << " "
      << std::setw(10) << a.get_y() << " "
      << std::setw(10) << a.get_z() << " "
      << std::setw(10) << a.get_a();

    return c;
}