#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>

namespace CXXCreator {

    gemmi::Structure CXXCreator_gemmi(const std::string& pdb_path) {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
        return st;
    }

} // namespace CXXCreator