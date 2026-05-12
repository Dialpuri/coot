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
};

}} // namespace coot::minimol