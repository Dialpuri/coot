#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <utility>
#include "coot-utils/atom-overlaps.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__link_atoms__96989c/gemmi/function.hh"

namespace coot {

inline std::string pad_atom_name_4(const std::string& name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    if (name.size() == 3) return " " + name;
    return name;
}

bool is_angle_related_via_link_gemmi(
    gemmi::CRA cra1,
    gemmi::CRA cra2,
    gemmi::Model* model_p,
    const std::vector<gemmi::Connection>& connections,
    const std::vector<std::pair<std::string, std::string>>& bonds_for_at_1,
    const std::vector<std::pair<std::string, std::string>>& bonds_for_at_2
);

} // namespace coot