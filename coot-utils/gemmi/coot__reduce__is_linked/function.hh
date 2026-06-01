#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/seqid.hpp>
#include <gemmi/util.hpp>
#include <string>

namespace coot {
namespace reduce {

/// Check whether the named atom of the given residue is part of any
/// connection (link) in the structure.
///
/// gemmi equivalent of coot::reduce::is_linked().
/// Pass the full structure so that structure.connections can be scanned.
inline bool is_linked_gemmi(
    const std::string &atom_name,
    gemmi::CRA cra,
    const gemmi::Structure &structure)
{
    bool status = false;

    if (cra.chain && cra.residue) {
        std::string chain_id  = cra.chain->name;
        int res_no            = cra.residue->seqid.num.value;
        char ins_code_char   = cra.residue->seqid.icode;

        // Normalise insertion code: gemmi uses ' ' for "none"
        std::string ins_code(1, ins_code_char);

        // gemmi trims atom names (e.g. " SG " -> "SG"), so we trim the input too
        std::string trimmed_atom = gemmi::trim_str(atom_name);

        for (const auto &conn : structure.connections) {
            const auto &addr1 = conn.partner1;
            const auto &addr2 = conn.partner2;

            // Check partner1
            if (chain_id == addr1.chain_name &&
                res_no   == addr1.res_id.seqid.num.value &&
                ins_code == std::string(1, addr1.res_id.seqid.icode) &&
                trimmed_atom == addr1.atom_name) {
                status = true;
                break;
            }

            // Check partner2
            if (chain_id == addr2.chain_name &&
                res_no   == addr2.res_id.seqid.num.value &&
                ins_code == std::string(1, addr2.res_id.seqid.icode) &&
                trimmed_atom == addr2.atom_name) {
                status = true;
                break;
            }
        }
    }

    return status;
}

} // namespace reduce
} // namespace coot