#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {

namespace {
// gemmi's padded_name() left-pads (PDB column 13-16), but MMDB also
// right-pads to 4 characters. This helper produces the same 4-char
// representation that MMDB GetAtomName() returned.
inline std::string atom_name_4(gemmi::Atom const& atom) {
    std::string s = atom.padded_name();
    s.resize(4, ' ');
    return s;
}
} // namespace

class beam_in_linked_residue {
public:
    // Ported from the MMDB version that takes mmdb::Residue* and returns
    // std::vector<mmdb::Atom*>.  The logic is identical: for each requested
    // atom name, scan the residue's atoms and push the first match.
    std::vector<gemmi::Atom*>
    get_atoms_gemmi(gemmi::Residue* residue,
                    const std::vector<std::string>& names) const {
        std::vector<gemmi::Atom*> v;
        for (const auto& name : names) {
            for (auto& atom : residue->atoms) {
                if (atom_name_4(atom) == name) {
                    v.push_back(&atom);
                    break;
                }
            }
        }
        return v;
    }
};

}  // namespace coot