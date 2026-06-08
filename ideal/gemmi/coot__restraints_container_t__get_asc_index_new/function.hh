#pragma once

#include <string>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__get_asc_index_old/gemmi/function.hh"

namespace coot {

inline int get_asc_index_new_gemmi(
    coot::restraints_container_t_gemmi& restraints,
    gemmi::Structure* st,
    const char* at_name,
    const char* alt_loc,
    int resno,
    const char* ins_code,
    const char* chain_id)
{
    (void)alt_loc;
    (void)ins_code;

    if (!st) {
        return -1;
    }

    // In gemmi there is no UDData mechanism, so the logic is:
    // 1. Find the atom in the structure (mirrors MMDB selection)
    // 2. If found, delegate to get_asc_index_old_gemmi
    // 3. If not found, return -1
    bool found = false;
    for (gemmi::Model& model : st->models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (gemmi::Residue& residue : chain.residues) {
                if (residue.seqid.num.value != resno) continue;
                for (gemmi::Atom& atom : residue.atoms) {
                    if (atom.name == at_name) {
                        found = true;
                        goto done_search;
                    }
                }
            }
        }
    }
done_search:
    if (!found) {
        return -1;
    }

    return restraints.get_asc_index_old_gemmi(at_name, resno, chain_id);
}

} // namespace coot
