#pragma once

#include <string>
#include <iostream>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__molecule__setup/gemmi/function.hh"

namespace coot {
namespace minimol {
class molecule;

inline int read_file_gemmi(molecule& mol, const std::string& pdb_filename) {
    try {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_filename);
        coot::minimol::setup_gemmi(mol, &st);
        return 0;
    } catch (const std::exception& e) {
        std::cout << "There was an error reading " << pdb_filename << ". \n";
        std::cout << "ERROR: " << e.what() << std::endl;
        return 15;
    }
}

} // namespace minimol
} // namespace coot