#include "function.hh"
#include <clipper/clipper.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace coot {
namespace smcif {

clipper::Cell get_cell_for_data_gemmi(const std::string& file_name) {
    clipper::Cell cell;
    
    std::ifstream file(file_name);
    if (!file.is_open()) {
        return cell;
    }
    
    std::string line;
    double a = 0.0, b = 0.0, c = 0.0;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    bool found_cell = false;
    
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Parse _cell_length_a, etc.
        if (line.find("_cell_length_a") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> a) found_cell = true;
        }
        else if (line.find("_cell_length_b") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> b) found_cell = true;
        }
        else if (line.find("_cell_length_c") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> c) found_cell = true;
        }
        else if (line.find("_cell_angle_alpha") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> alpha) found_cell = true;
        }
        else if (line.find("_cell_angle_beta") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> beta) found_cell = true;
        }
        else if (line.find("_cell_angle_gamma") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag;
            if (iss >> gamma) found_cell = true;
        }
    }
    
    file.close();
    
    if (found_cell && a > 0 && b > 0 && c > 0 && alpha > 0 && beta > 0 && gamma > 0) {
        clipper::Cell_descr cell_descr(a, b, c,
                                      clipper::Util::d2rad(alpha),
                                      clipper::Util::d2rad(beta),
                                      clipper::Util::d2rad(gamma));
        cell = clipper::Cell(cell_descr);
    }
    
    return cell;
}

} // namespace smcif
} // namespace coot