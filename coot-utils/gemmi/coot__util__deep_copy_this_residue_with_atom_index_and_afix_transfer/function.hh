#pragma once

#include <gemmi/model.hpp>
#include <memory>
#include <string>

namespace coot {
namespace util {

/// gemmi port of deep_copy_this_residue_with_atom_index_and_afix_transfer
///
/// Deep-copies a residue from the given CRA, filtering atoms by altloc.
/// UDD data transfer from the MMDB version is not applicable in gemmi
/// (gemmi has no UDD system).
///
/// @param cra                 Parent context (chain*, residue*, atom*).
/// @param altconf             Alternate conformation filter (ignored when
///                            whole_residue_flag is true).
/// @param whole_residue_flag  If true, copy all atoms regardless of altloc.
/// @return Newly-allocated residue, or nullptr if cra.residue is null.
std::unique_ptr<gemmi::Residue>
deep_copy_this_residue_with_atom_index_and_afix_transfer_gemmi(
    gemmi::CRA cra,
    const std::string &altconf,
    short int whole_residue_flag)
{
    if (!cra.residue) {
        return nullptr;
    }

    auto new_res = std::make_unique<gemmi::Residue>();
    new_res->name  = cra.residue->name;
    new_res->seqid = cra.residue->seqid;

    for (const auto& atom : cra.residue->atoms) {
        bool copy = false;
        if (whole_residue_flag) {
            copy = true;
        } else {
            // In gemmi: '\0' = no alternate conformation (MMDB "")
            // Match if atom has no altloc, or altloc equals the first
            // character of altconf.
            copy = (atom.altloc == '\0' ||
                    (!altconf.empty() && atom.altloc == altconf[0]));
        }
        if (copy) {
            new_res->atoms.push_back(atom);
        }
    }

    return new_res;
}

} // namespace util
} // namespace coot