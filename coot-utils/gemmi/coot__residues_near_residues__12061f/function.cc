#include "function.hh"
#include <gemmi/model.hpp>
#include <vector>

namespace coot {

namespace {

std::string icode_to_string(char icode) {
    if (icode != ' ' && icode != '\0')
        return std::string(1, icode);
    return std::string();
}

residue_spec_t spec_from_gemmi(const gemmi::Chain& chain, const gemmi::Residue& res, int model_number) {
    return residue_spec_t(model_number, chain.name, res.seqid.num.value, icode_to_string(res.seqid.icode));
}

} // anon namespace

std::map<coot::residue_spec_t, std::set<coot::residue_spec_t>>
residues_near_residues_gemmi(gemmi::Structure* st, float dist_crit) {
    std::map<coot::residue_spec_t, std::set<coot::residue_spec_t>> result;

    if (st == nullptr) {
        return result;
    }

    if (st->models.empty()) {
        return result;
    }

    if (dist_crit <= 0.0f) {
        return result;
    }

    gemmi::Model& model = st->models[0];

    // Collect all atoms with their chain/residue pointers (mimicking atom_selection)
    struct AtomEntry {
        const gemmi::Atom* atom;
        const gemmi::Chain* chain;
        const gemmi::Residue* residue;
    };
    std::vector<AtomEntry> atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (const gemmi::Atom& atom : res.atoms) {
                atoms.push_back({&atom, &chain, &res});
            }
        }
    }

    double radius_sq = static_cast<double>(dist_crit) * static_cast<double>(dist_crit);

    // For each pair (i, j), check distance and add mapping
    // The original MMDB SeekContacts returns contacts bidirectionally
    // (both A->B and B->A for each pair), so we iterate all ordered pairs
    // where i != j to match that behavior.
    size_t n = atoms.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;

            double dx = atoms[i].atom->pos.x - atoms[j].atom->pos.x;
            double dy = atoms[i].atom->pos.y - atoms[j].atom->pos.y;
            double dz = atoms[i].atom->pos.z - atoms[j].atom->pos.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;

            if (dist_sq <= radius_sq) {
                auto key_spec = spec_from_gemmi(*atoms[i].chain, *atoms[i].residue, 0);
                auto data_spec = spec_from_gemmi(*atoms[j].chain, *atoms[j].residue, 0);

                // Only add if residues are different (matching original logic: if (data != key))
                if (data_spec != key_spec) {
                    result[key_spec].insert(data_spec);
                }
            }
        }
    }

    return result;
}

} // namespace coot