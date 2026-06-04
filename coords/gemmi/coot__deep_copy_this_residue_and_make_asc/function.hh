#pragma once

#include <gemmi/model.hpp>
#include <memory>
#include <string>
#include "coot-utils/gemmi/make_asc/function.hh"

namespace coot {

std::pair<std::unique_ptr<gemmi::Residue>, atom_selection_container_t>
deep_copy_this_residue_and_make_asc_gemmi(
    const gemmi::Structure& orig_st,
    gemmi::CRA cra,
    const std::string& altconf,
    short int whole_residue_flag,
    int atom_index_handle,
    int udd_afix_handle);

}  // namespace coot