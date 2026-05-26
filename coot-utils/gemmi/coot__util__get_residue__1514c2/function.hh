#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_by_binary_search/gemmi/function.hh"

namespace coot {
namespace util {

// Gemmi version of get_residue — returns residue via CRA (Chain-Residue-Atom)
// If not found, returns {nullptr, nullptr, nullptr}
inline gemmi::CRA get_residue_gemmi(const std::string &chain_id,
                                    int reso, const std::string &insertion_code,
                                    gemmi::Structure &st) {
    if (st.models.empty())
        return {nullptr, nullptr, nullptr};

    // Direct linear search over chains and residues (mimics original MMDB logic)
    std::string icode_str = insertion_code.empty() ? " " : insertion_code;
    char ins_code_char = icode_str.empty() ? ' ' : icode_str[0];

    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto &residue : chain.residues) {
                    if (residue.seqid.num.value == reso) {
                        if (residue.seqid.icode == ins_code_char) {
                            return {&chain, &residue, nullptr};
                        }
                    }
                }
            }
        }
    }
    return {nullptr, nullptr, nullptr};
}

} // namespace util
} // namespace coot