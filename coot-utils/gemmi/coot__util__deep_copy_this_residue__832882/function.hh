#pragma once

#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

inline gemmi::Residue *
deep_copy_this_residue_gemmi(const gemmi::Residue *residue,
                             const std::pair<bool, std::string> &use_alt_conf) {
    gemmi::Residue *rres = nullptr;

    if (residue) {
        rres = new gemmi::Residue;
        rres->seqid.num.value = residue->seqid.num.value;
        rres->name = residue->name;
        rres->seqid.icode = residue->seqid.icode;

        for (const gemmi::Atom& atom : residue->atoms) {
            if (use_alt_conf.first) {
                // In MMDB, altLoc is a char array. Empty means no altloc (blank).
                // In gemmi, altloc is a char. If it's ' ' (space) or '\0', treat as empty.
                if (atom.altloc != ' ' && atom.altloc != '\0') {
                    if (std::string(1, atom.altloc) != use_alt_conf.second)
                        continue;
                }
            }
            gemmi::Atom atom_p = atom;
            rres->atoms.push_back(atom_p);
        }
    }

    return rres;
}

} // namespace util
} // namespace coot