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

    atom() {
        int_user_data = -1;
        occupancy = -1;
        temperature_factor = -1;
        name = "";
        element = "";
        altLoc = "";
        pos = gemmi::Vec3(0.0, 0.0, 0.0);
    }
};

// gemmi port of the MMDB constructor:
// atom::atom(std::string atom_name, std::string ele, const clipper::Coord_orth &pos_in,
//            const std::string &altloc, float dbf)
// clipper::Coord_orth → gemmi::Vec3
atom atom_gemmi(std::string atom_name,
                std::string ele,
                gemmi::Vec3 pos_in,
                const std::string &altloc,
                float dbf);

}} // namespace coot::minimol
