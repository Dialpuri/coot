#pragma once

#include <string>
#include <gemmi/math.hpp>

namespace coot { namespace minimol {

struct atom {
    std::string name;
    std::string element;
    std::string altLoc;
    gemmi::Vec3 pos;
    float occupancy;
    float temperature_factor;
    int int_user_data;

    // Default constructor
    atom()
        : pos(0.0, 0.0, 0.0),
          occupancy(-1),
          temperature_factor(-1),
          int_user_data(-1) {}

    // gemmi port of the MMDB constructor:
    // atom::atom(std::string atom_name, std::string ele,
    //            const clipper::Coord_orth &pos_in, const std::string &altloc,
    //            float occupancy_in, float dbf)
    static atom atom_gemmi(std::string atom_name,
                           std::string ele,
                           gemmi::Vec3 pos_in,
                           const std::string &altloc,
                           float occupancy_in,
                           float dbf) {
        atom result;
        result.name  = atom_name;
        result.element = ele;
        result.pos   = pos_in;
        result.altLoc = altloc;
        result.occupancy = occupancy_in;
        result.temperature_factor = dbf;
        result.int_user_data = -1;
        return result;
    }
};

}} // namespace coot::minimol
