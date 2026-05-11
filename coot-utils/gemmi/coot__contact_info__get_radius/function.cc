#include "function.hh"
#include <algorithm>

namespace coot {

void contact_info::setup_atom_radii() {
    atom_radii.resize(23);
    atom_radii[0]  = std::make_pair(" C", 0.77);
    atom_radii[1]  = std::make_pair(" N", 0.65);
    atom_radii[2]  = std::make_pair(" O", 0.6);
    atom_radii[3]  = std::make_pair(" H", 0.35);
    atom_radii[4]  = std::make_pair(" S", 1.1);
    atom_radii[5]  = std::make_pair(" P", 1.0);
    atom_radii[6]  = std::make_pair("SE", 1.15);
    atom_radii[7]  = std::make_pair("BR", 1.15);
    atom_radii[8]  = std::make_pair("CL", 1.0);
    atom_radii[9]  = std::make_pair(" I", 1.4);
    atom_radii[10] = std::make_pair(" F", 0.5);
    atom_radii[11] = std::make_pair(" K", 2.2);
    atom_radii[12] = std::make_pair("AS", 1.3);
    atom_radii[13] = std::make_pair("NA", 1.8);
    atom_radii[14] = std::make_pair("MG", 1.5);
    atom_radii[15] = std::make_pair("AU", 1.4);
    atom_radii[16] = std::make_pair("BE", 1.05);
    atom_radii[17] = std::make_pair("FE", 1.4);
    atom_radii[18] = std::make_pair("ZN", 1.35);
    atom_radii[19] = std::make_pair("PD", 1.6);
    atom_radii[20] = std::make_pair("PB", 1.46);
    atom_radii[21] = std::make_pair("PT", 1.46);
    atom_radii[22] = std::make_pair("AG", 1.36);
}

double contact_info::get_radius_gemmi(const std::string &element) const {
    double r = 0.9;
    for (size_t i = 0; i < atom_radii.size(); ++i) {
        if (atom_radii[i].first == element) {
            r = atom_radii[i].second;
            break;
        }
    }
    return r;
}

} // namespace coot