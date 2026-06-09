#include "function.hh"

namespace coot {
namespace h_bonds {

std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>>
make_neighbour_map_gemmi(const std::vector<gemmi::CRA>& cras_1,
                         const std::vector<gemmi::CRA>& cras_2) {
    std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>> atom_map;

    // Helper: squared distance between two atoms
    auto sq_dist = [](gemmi::Atom const& a, gemmi::Atom const& b) {
        double dx = a.pos.x - b.pos.x;
        double dy = a.pos.y - b.pos.y;
        double dz = a.pos.z - b.pos.z;
        return static_cast<float>(std::sqrt(dx*dx + dy*dy + dz*dz));
    };

    // First pass: atoms from cras_1, distance 0.1–1.8 Å
    for (size_t i = 0; i < cras_1.size(); ++i) {
        if (!cras_1[i].atom) continue;
        for (size_t j = 0; j < cras_1.size(); ++j) {
            if (i == j) continue;
            if (!cras_1[j].atom) continue;
            // Both must be in the same residue
            if (cras_1[i].residue != cras_1[j].residue) continue;
            float d = sq_dist(*cras_1[i].atom, *cras_1[j].atom);
            if (d >= 0.1f && d <= 1.8f) {
                atom_map[cras_1[i].atom].push_back({cras_1[j].atom, d});
            }
        }
    }

    // Second pass: atoms from cras_2, distance 0.1–1.7 Å, deduped
    for (size_t i = 0; i < cras_2.size(); ++i) {
        if (!cras_2[i].atom) continue;
        for (size_t j = 0; j < cras_2.size(); ++j) {
            if (i == j) continue;
            if (!cras_2[j].atom) continue;
            // Both must be in the same residue
            if (cras_2[i].residue != cras_2[j].residue) continue;
            float d = sq_dist(*cras_2[i].atom, *cras_2[j].atom);
            if (d >= 0.1f && d <= 1.7f) {
                // Only add if not already present
                auto& vec = atom_map[cras_2[i].atom];
                bool found = false;
                for (const auto& p : vec) {
                    if (p.first == cras_2[j].atom) { found = true; break; }
                }
                if (!found) {
                    vec.push_back({cras_2[j].atom, d});
                }
            }
        }
    }

    // Sort neighbours by distance (closest first)
    for (auto& entry : atom_map) {
        std::sort(entry.second.begin(), entry.second.end(),
                  [](const std::pair<gemmi::Atom*, float>& a,
                     const std::pair<gemmi::Atom*, float>& b) {
                      return a.second < b.second;
                  });
    }

    return atom_map;
}

} // namespace h_bonds
} // namespace coot
