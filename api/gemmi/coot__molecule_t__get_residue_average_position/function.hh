#pragma once
#include <string>
#include <vector>
#include <utility>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Parse CID format "//chain/resnum" and return {chain_name, residue_number}
// Returns empty string for chain if parsing fails
inline std::pair<std::string, int> parse_cid(const std::string& cid) {
    // CID format is "//chain/resnum"
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/') {
        return {"", 0};
    }
    
    size_t slash_pos = cid.find('/', 2);
    if (slash_pos == std::string::npos) {
        return {"", 0};
    }
    
    std::string chain = cid.substr(2, slash_pos - 2);
    int resnum = 0;
    try {
        resnum = std::stoi(cid.substr(slash_pos + 1));
    } catch (...) {
        return {"", 0};
    }
    
    return {chain, resnum};
}

// Get average position of all atoms in a residue
// Returns empty vector if residue not found
// Returns {x, y, z} coordinates of average position
std::vector<double> get_residue_average_position_gemmi(
    const gemmi::Structure& st,
    const std::string& cid
);

} // namespace molecule_t
} // namespace coot