#pragma once
#include <gemmi/model.hpp>
#include <utility>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__max_resno_in_chain/gemmi/function.hh"

namespace coot {
namespace util {

std::pair<bool, int> max_resno_in_molecule_gemmi(const gemmi::Structure& st) {
    bool found_residues = false;
    int current_high = -31999;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            std::pair<bool, int> p = max_resno_in_chain_gemmi(chain);
            if (p.first) {
                if (p.second > current_high) {
                    current_high = p.second;
                    found_residues = true;
                }
            }
        }
    }
    return std::pair<bool, int>(found_residues, current_high);
}

} // namespace util
} // namespace coot