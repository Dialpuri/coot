#pragma once

#include <vector>
#include <string>
#include <utility>

namespace coot {

std::vector<std::pair<std::string, double> > setup_atom_radii_gemmi() {
    std::vector<std::pair<std::string, double> > atom_radii;
    atom_radii.resize(23);
    atom_radii[ 0] = std::pair<std::string, double> (" C", 0.77);
    atom_radii[ 1] = std::pair<std::string, double> (" N", 0.65);
    atom_radii[ 2] = std::pair<std::string, double> (" O", 0.6);
    atom_radii[ 3] = std::pair<std::string, double> (" H", 0.35);
    // atom_radii[ 4] = std::pair<std::string, double> (" S", 0.9);
    atom_radii[ 4] = std::pair<std::string, double> (" S", 1.1); // S-S bonds 2.16A?
    atom_radii[ 5] = std::pair<std::string, double> (" P", 1.0);
    atom_radii[ 6] = std::pair<std::string, double> ("SE", 1.15);
    atom_radii[ 7] = std::pair<std::string, double> ("BR", 1.15);
    atom_radii[ 8] = std::pair<std::string, double> ("CL", 1.0);
    atom_radii[ 9] = std::pair<std::string, double> (" I", 1.4);
    atom_radii[10] = std::pair<std::string, double> (" F", 0.5);
    atom_radii[11] = std::pair<std::string, double> (" K", 2.2);
    atom_radii[12] = std::pair<std::string, double> ("AS", 1.3);
    atom_radii[13] = std::pair<std::string, double> ("NA", 1.8);
    atom_radii[14] = std::pair<std::string, double> ("MG", 1.5);
    atom_radii[15] = std::pair<std::string, double> ("AU", 1.4);
    atom_radii[16] = std::pair<std::string, double> ("BE", 1.05);
    atom_radii[17] = std::pair<std::string, double> ("FE", 1.4);
    atom_radii[18] = std::pair<std::string, double> ("ZN", 1.35);
    atom_radii[19] = std::pair<std::string, double> ("PD", 1.6);
    atom_radii[20] = std::pair<std::string, double> ("PB", 1.46);
    atom_radii[21] = std::pair<std::string, double> ("PT", 1.46);
    atom_radii[22] = std::pair<std::string, double> ("AG", 1.36);
    return atom_radii;
}

}