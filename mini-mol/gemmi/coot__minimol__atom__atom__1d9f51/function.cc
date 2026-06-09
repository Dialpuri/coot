#include "function.hh"

namespace coot { namespace minimol {

atom atom_gemmi(std::string atom_name,
                std::string ele,
                gemmi::Vec3 pos_in,
                const std::string &altloc,
                float dbf) {
    atom at;
    at.name = atom_name;
    at.element = ele;
    at.pos = pos_in;
    at.altLoc = altloc;
    at.occupancy = 1.0f;
    at.temperature_factor = dbf;
    at.int_user_data = -1;
    return at;
}

}} // namespace coot::minimol
