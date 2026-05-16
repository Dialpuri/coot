#pragma once

#include <map>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace side_chain_densities {

// Ported from coot::side_chain_densities::make_sequence_for_chain (MMDB).
// Returns a map from residue sequence number to 3-letter residue name
// for every residue in the given gemmi chain.
std::map<int, std::string>
make_sequence_for_chain_gemmi(const gemmi::Chain &chain);

} // namespace side_chain_densities
} // namespace coot

inline
std::map<int, std::string>
coot::side_chain_densities::make_sequence_for_chain_gemmi(const gemmi::Chain &chain) {
    std::map<int, std::string> m;
    for (const gemmi::Residue &residue : chain.residues) {
        int res_no = residue.seqid.num.value;
        std::string res_name = residue.name;
        m[res_no] = res_name;
    }
    return m;
}