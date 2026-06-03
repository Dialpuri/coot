#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>

namespace coot {

class ShelxIns {
public:
    std::map<std::string, unsigned int> get_atomic_contents_gemmi(const gemmi::Structure& st) const {
        std::map<std::string, unsigned int> m;
        if (st.models.empty()) return m;
        const gemmi::Model& model = st.models[0];
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                for (const gemmi::Atom& atom : residue.atoms) {
                    std::string ele = atom.element.name();
                    // MMDB stores raw PDB element string (right-justified in 2 chars):
                    // single-char elements get a leading space: " C", " N", etc.
                    if (ele.size() == 1) {
                        ele = " " + ele;
                    }
                    if (!ele.empty()) {
                        m[ele]++;
                    }
                }
            }
        }
        return m;
    }
};

} // namespace coot