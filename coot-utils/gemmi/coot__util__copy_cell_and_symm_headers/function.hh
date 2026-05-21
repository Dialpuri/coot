#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

bool copy_cell_and_symm_headers_gemmi(gemmi::Structure& src, gemmi::Structure& dest) {
    if (&src && &dest) {
        dest.cell = src.cell;
        dest.spacegroup_hm = src.spacegroup_hm;
        return true;
    }
    return false;
}

}}