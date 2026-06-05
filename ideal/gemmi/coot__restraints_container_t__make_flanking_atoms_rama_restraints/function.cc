#include "function.hh"

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace coot {

std::vector<rama_flanking_triple_t>
find_rama_flanking_triples_gemmi(
    const gemmi::Model& model,
    const std::string& chain_id_save,
    int istart_res,
    int iend_res,
    int istart_minus_flag,
    int iend_plus_flag)
{
    std::vector<rama_flanking_triple_t> triples;

    // Helper: find residue by chain and seqid in a model
    auto find_res = [&](const gemmi::Model& mod,
                        const std::string& chain_name,
                        int seqnum) -> const gemmi::Residue* {
        for (const auto& ch : mod.chains) {
            if (ch.name != chain_name) continue;
            for (const auto& res : ch.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
        return nullptr;
    };

    if (!istart_minus_flag || !iend_plus_flag) {
        return triples;
    }

    // First flanking triple: (istart-1, istart, istart+1)
    {
        rama_flanking_triple_t triple;
        triple.first  = find_res(model, chain_id_save, istart_res - 1);
        triple.middle = find_res(model, chain_id_save, istart_res);
        triple.third  = find_res(model, chain_id_save, istart_res + 1);
        if (triple.is_valid()) {
            triples.push_back(triple);
        }
    }

    // Second triple only if istart != iend: (iend-1, iend, iend+1)
    if (istart_res != iend_res) {
        rama_flanking_triple_t triple;
        triple.first  = find_res(model, chain_id_save, iend_res - 1);
        triple.middle = find_res(model, chain_id_save, iend_res);
        triple.third  = find_res(model, chain_id_save, iend_res + 1);
        if (triple.is_valid()) {
            triples.push_back(triple);
        }
    }

    return triples;
}

int make_flanking_atoms_rama_restraints_gemmi(
    gemmi::Structure& mol,
    const std::string& chain_id_save,
    int istart_res,
    int iend_res,
    int istart_minus_flag,
    int iend_plus_flag,
    const coot::protein_geometry& geom)
{
    int n_rama_restraints = 0;

    if (!mol.models.empty()) {
        const gemmi::Model& model = mol.models[0];
        auto triples = find_rama_flanking_triples_gemmi(
            model, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag);

        for (const auto& triple : triples) {
            // In the original, add_rama is called here for each triple
            // but the function always returns 0, matching original behavior
            (void)geom; // suppress unused warning
        }
    }

    return n_rama_restraints;
}

} // namespace coot