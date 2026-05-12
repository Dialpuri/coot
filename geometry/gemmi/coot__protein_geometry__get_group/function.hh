#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

struct protein_geometry {
    std::string get_group_gemmi(const gemmi::Residue& residue) const;
    std::string get_group(const std::string& res_name) const;
};

} // namespace coot