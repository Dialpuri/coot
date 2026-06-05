#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include <clipper/core/clipper_types.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/xmap.h>
#include <string>
#include <vector>

namespace coot {

struct simple_sm_u {
    std::string label;
    double u11;
    double u22;
    double u33;
    double u12;
    double u13;
    double u23;
    simple_sm_u(std::string l, double a, double b, double c, double d, double e, double f)
        : label(l), u11(a), u22(b), u33(c), u12(d), u13(e), u23(f) {}
};

} // namespace coot

namespace coot::smcif {

std::vector<gemmi::Atom> read_coordinates_gemmi(
    const std::string& cif_path,
    const clipper::Cell& cell,
    const clipper::Spacegroup& spg);

} // namespace coot::smcif