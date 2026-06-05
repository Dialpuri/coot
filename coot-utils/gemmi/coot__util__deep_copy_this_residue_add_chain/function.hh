#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>

namespace coot {
namespace util {

/*
 * gemmi port of coot::util::deep_copy_this_residue_add_chain
 *
 * Deep-copies a residue (filtered by alt-loc conformer), optionally attaching
 * it to a newly-allocated chain.
 *
 * Returns a pair<Residue*, Chain*>:
 *   - residue pointer is always valid for non-null input
 *   - chain pointer is non-null only when attach_to_new_chain_flag is true
 *     (and caller owns the chain, which in turn owns the residue)
 *   - when attach_to_new_chain_flag is false, the caller owns the residue
 */
inline std::pair<gemmi::Residue*, gemmi::Chain*> deep_copy_this_residue_add_chain_gemmi(
    gemmi::CRA cra,
    const std::string& altconf,
    bool whole_residue_flag,
    bool attach_to_new_chain_flag)
{
    gemmi::Residue* rres = nullptr;
    gemmi::Chain* chain_p = nullptr;

    if (cra.residue) {
        if (attach_to_new_chain_flag) {
            chain_p = new gemmi::Chain(cra.chain ? cra.chain->name : "");
        }

        rres = new gemmi::Residue;
        rres->name = cra.residue->name;
        rres->seqid.num.value = cra.residue->seqid.num.value;
        rres->seqid.icode = cra.residue->seqid.icode;

        // altconf is a string; in gemmi altloc is a char ('\0' means "no alt loc")
        char altconf_char = altconf.empty() ? '\0' : altconf[0];

        for (const auto& atom : cra.residue->atoms) {
            // gemmi has no terminator-atom concept, so the original isTer() guard
            // is a no-op — all atoms pass.
            if (whole_residue_flag ||
                atom.altloc == altconf_char ||
                atom.altloc == '\0') {
                rres->atoms.push_back(atom);
            }
        }

        if (attach_to_new_chain_flag) {
            chain_p->residues.push_back(*rres);
            delete rres;
            rres = &chain_p->residues.back();
        }
    }

    return {rres, chain_p};
}

} // namespace util
} // namespace coot