#pragma once
#include <string>
#include <vector>
#include <map>
#include <gemmi/model.hpp>

// Include the 7-param generated port — our 6-param wrapper delegates to it
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_amino_hydrogens__2025be/gemmi/function.hh"

namespace coot {
namespace reduce {

/// gemmi port of coot::reduce::add_amino_hydrogens (6-param / map-based version).
///
/// Resolves the first second-neighbor and the first third-neighbor from the
/// supplied vectors & maps, then delegates to the 7-param overload.
///
/// `structure` and `chain` are forwarded to the 7-param delegate for neighbor
/// lookups; pass nullptr for chain if the caller does not have chain context.
inline void add_amino_hydrogens_gemmi(
    const std::string &H_at_name_1,
    const std::string &H_at_name_2,
    const std::string &first_neighb,
    const std::vector<std::string> &second_neighb_vec,
    const std::map<std::string, std::vector<std::string>> &third_neighb_map,
    double bl_amino,
    gemmi::Residue &residue,
    const gemmi::Structure &structure,
    gemmi::Chain *chain = nullptr)
{
    if (!second_neighb_vec.empty()) {
        std::string second = second_neighb_vec[0];
        auto it = third_neighb_map.find(second);
        if (it != third_neighb_map.end()) {
            const auto& thirds = it->second;
            if (!thirds.empty()) {
                const std::string& third = thirds[0];
                // Dispatch to the 7-param overload (from included header)
                coot::reduce::add_amino_hydrogens_gemmi(H_at_name_1, H_at_name_2,
                                          first_neighb, second, third,
                                          bl_amino, residue, structure, chain);
            }
        }
    }
}

} // namespace reduce
} // namespace coot