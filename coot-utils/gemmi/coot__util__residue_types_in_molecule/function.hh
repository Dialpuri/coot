#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

inline std::vector<std::string> residue_types_in_molecule_gemmi(const gemmi::Structure& st) {
    std::vector<std::string> v;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                std::string resname = residue.name;
                bool found = false;
                for (const auto& existing : v) {
                    if (existing == resname) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    v.push_back(resname);
                }
            }
        }
    }
    
    return v;
}

} // namespace util
} // namespace coot