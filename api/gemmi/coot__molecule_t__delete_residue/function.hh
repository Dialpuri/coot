#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__pdbcleanup_serial_residue_numbers/gemmi/function.hh"

namespace coot {

/// Delete a residue identified by residue_spec from the given model.
/// Returns 1 if a residue was deleted, 0 otherwise.
/// Cleans up residue numbering after deletion.
inline int delete_residue_gemmi(gemmi::Structure& st, gemmi::Model& model, coot::residue_spec_t const& residue_spec) {
    int was_deleted = 0;

    // Normalize insertion code: MMDB uses "" for "no insertion code", gemmi uses ' '
    auto norm = [](const std::string& ic) { return ic.empty() ? std::string(" ") : ic; };

    for (auto& chain : model.chains) {
        if (residue_spec.chain_id == chain.name) {
            for (auto it = chain.residues.begin(); it != chain.residues.end(); ++it) {
                if (it->seqid.num.value == residue_spec.res_no) {
                    std::string gemmi_icode(1, it->seqid.icode);
                    if (norm(residue_spec.ins_code) == gemmi_icode) {
                        it = chain.residues.erase(it);
                        was_deleted = 1;
                        break;
                    }
                }
            }
            if (was_deleted) break;
        }
    }

    if (was_deleted) {
        coot::util::pdbcleanup_serial_residue_numbers_gemmi(st);
    }

    return was_deleted;
}

} // namespace coot