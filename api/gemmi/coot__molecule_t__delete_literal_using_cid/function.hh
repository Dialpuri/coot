#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/polyheur.hpp>

namespace coot {

struct molecule_t {
    gemmi::Structure structure;
    
    molecule_t() = default;
    explicit molecule_t(gemmi::Structure str) : structure(std::move(str)) {}
    
    int delete_literal_using_cid_gemmi(const std::string &atom_selection_cids);
};

} // namespace coot