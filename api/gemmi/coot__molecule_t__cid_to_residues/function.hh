#pragma once

#include <string>
#include <vector>
#include <set>
#include <gemmi/model.hpp>

namespace coot {

inline std::vector<std::string> split_string(const std::string &s, const std::string &delim) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = s.find(delim, start);
    while (end != std::string::npos) {
        result.push_back(s.substr(start, end - start));
        start = end + delim.size();
        end = s.find(delim, start);
    }
    result.push_back(s.substr(start));
    return result;
}

std::vector<gemmi::Residue *>
cid_to_residues_gemmi(const std::string &atom_selection_cids, gemmi::Structure &st);

}  // namespace coot