#pragma once
#include <string>
#include <gemmi/math.hpp>

namespace coot { namespace minimol {

struct atom {
    std::string name;
    std::string element;
    std::string altLoc;
    gemmi::Position pos;
    float occupancy;
    float temperature_factor;
    int int_user_data;

    atom() {
        int_user_data = -1;
        occupancy = -1;
        temperature_factor = -1;
        name = "";
        element = "";
        altLoc = "";
        pos = gemmi::Position(0.0, 0.0, 0.0);
    }

    // gemmi port of the MMDB constructor:
    // atom::atom(std::string atom_name, std::string ele, float x, float y, float z,
    //            const std::string &altloc, float occupancy_in, float dbf)
    static atom atom_gemmi(std::string atom_name,
                           std::string ele,
                           float x, float y, float z,
                           const std::string &altloc,
                           float occupancy_in,
                           float dbf) {
        atom result;
        result.name = atom_name;
        result.element = ele;
        result.altLoc = altloc;
        result.pos = gemmi::Position(x, y, z);
        result.occupancy = occupancy_in;
        result.temperature_factor = dbf;
        result.int_user_data = -1;
        return result;
    }
};

}} // namespace coot::minimol
