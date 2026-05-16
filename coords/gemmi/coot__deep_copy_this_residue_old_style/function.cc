#include "function.hh"
#include <gemmi/model.hpp>
#include <string>

namespace coot {

gemmi::Residue *
deep_copy_this_residue_old_style_gemmi(const gemmi::CRA &cra,
                                       const std::string &altconf,
                                       int whole_residue_flag,
                                       bool embed_in_chain_flag) {

    const gemmi::Residue &res = *cra.residue;
    gemmi::Residue *rres = nullptr;

    if (!res.atoms.empty()) {
        rres = new gemmi::Residue();
        rres->name = res.name;
        rres->seqid = res.seqid;

        int n_added_atoms = 0;
        for (const auto &atom : res.atoms) {
            // Skip TER atoms (equivalent to at->isTer() in MMDB)
            if (atom.name == "TER") {
                continue;
            }
            // In MMDB: altLoc is a string. "" means no altloc.
            // In gemmi: altloc is a char. ' ' means no altloc.
            bool is_no_altloc = (atom.altloc == ' ' || atom.altloc == '\0');

            if (whole_residue_flag) {
                gemmi::Atom new_atom = atom;
                rres->atoms.push_back(new_atom);
                n_added_atoms++;
            } else if (altconf.empty()) {
                // Match atoms with no altloc (altLoc == "" in MMDB)
                if (is_no_altloc) {
                    gemmi::Atom new_atom = atom;
                    rres->atoms.push_back(new_atom);
                    n_added_atoms++;
                }
            } else {
                // Match atoms with the specific altloc character
                if (atom.altloc == altconf[0]) {
                    gemmi::Atom new_atom = atom;
                    rres->atoms.push_back(new_atom);
                    n_added_atoms++;
                }
            }
        }

        if (n_added_atoms == 0) {
            delete rres;
            rres = nullptr;
        }
    }

    return rres;
}

} // namespace coot