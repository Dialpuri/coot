#pragma once

#include <gemmi/model.hpp>
#include <set>
#include <string>
#include <algorithm>

namespace coot::util {

bool delete_all_carbohydrate_gemmi(gemmi::Structure* st) {
    static const std::set<std::string> cho_set = {
        "NAG", "MAN", "BMA", "FUL",
        "FUC", "XYP", "SIA", "GAL",
        "NDG", "BGC", "A2G"
    };

    if (!st) return false;

    bool deleted = false;

    for (gemmi::Model& model : st->models) {
        for (gemmi::Chain& chain : model.chains) {
            auto it = std::remove_if(
                chain.residues.begin(),
                chain.residues.end(),
                [&](const gemmi::Residue& res) {
                    return cho_set.find(res.name) != cho_set.end();
                }
            );
            std::size_t removed_count = std::distance(it, chain.residues.end());
            chain.residues.erase(it, chain.residues.end());
            if (removed_count > 0) deleted = true;
        }
    }

    return deleted;
}

} // namespace coot::util