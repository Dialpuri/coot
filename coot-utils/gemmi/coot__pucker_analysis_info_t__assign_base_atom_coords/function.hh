#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

struct pucker_analysis_info_t {
    std::vector<gemmi::Vec3> base_atoms_coords;
    gemmi::Atom* N1_or_9 = nullptr;
    gemmi::Atom* C1_prime = nullptr;
    std::string altconf;

    pucker_analysis_info_t(gemmi::Residue* res, const std::string& alt)
        : altconf(alt) {
        assign_base_atom_coords_gemmi(res);
        if (base_atoms_coords.size() < 3) {
            throw std::runtime_error("base lsq plane has no value");
        }
    }

    void assign_base_atom_coords_gemmi(gemmi::Residue* residue);
};

} // namespace coot