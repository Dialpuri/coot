#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {
namespace util {

class backrub_residue_triple_t {
public:
    gemmi::Residue* this_residue;
    gemmi::Residue* prev_residue;
    gemmi::Residue* next_residue;
    std::string alt_conf;

    backrub_residue_triple_t(gemmi::Residue* prev_r, gemmi::Residue* this_r,
                             gemmi::Residue* next_r, std::string alt)
        : prev_residue(prev_r), this_residue(this_r),
          next_residue(next_r), alt_conf(std::move(alt)) {}

    void trim_residue_atoms_generic_gemmi(gemmi::Residue* residue_p,
                                          std::vector<std::string> keep_atom_vector,
                                          bool use_keep_atom_vector);
};

} // namespace util
} // namespace coot
