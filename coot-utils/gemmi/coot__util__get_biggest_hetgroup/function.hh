#pragma once

#include <gemmi/model.hpp>
#include <gemmi/to_pdb.hpp>
#include <vector>

namespace coot { namespace util {

/// gemmi port of coot::util::get_biggest_hetgroup.
/// Returns the CRA (Chain*, Residue*, nullptr) of the hetero residue with
/// the most atoms in the given model (matching MMDB GetModel(1)).
/// Returns CRA{nullptr,nullptr,nullptr} when model is nullptr or there are
/// no hetero residues.
gemmi::CRA get_biggest_hetgroup_gemmi(gemmi::Model* model);

inline gemmi::CRA get_biggest_hetgroup_gemmi(gemmi::Model* model) {
    if (!model) return {nullptr, nullptr, nullptr};

    gemmi::CRA best{};
    int biggest_n = 0;

    for (auto& chain : model->chains) {
        for (auto& residue : chain.residues) {
            if (gemmi::use_hetatm(residue)) {
                int n_atoms = static_cast<int>(residue.atoms.size());
                if (n_atoms > biggest_n) {
                    best = {&chain, &residue, nullptr};
                    biggest_n = n_atoms;
                }
            }
        }
    }

    return best;
}

}} // namespace coot::util