#include "function.hh"
#include "coot/coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

static bool is_water_or_hoh(const gemmi::Residue &res) {
    std::string resname = res.name;
    return (resname == "WAT" || resname == "HOH");
}

int trim_molecule_by_map_gemmi(gemmi::Model &model,
                               const clipper::Xmap<float> &xmap,
                               float map_level,
                               short int remove_or_zero_occ_flag,
                               short int waters_only_flag) {
    int n_changed = 0;

    for (gemmi::Chain &chain : model.chains) {
        for (gemmi::Residue &residue : chain.residues) {
            bool is_water = is_water_or_hoh(residue);

            if (waters_only_flag) {
                if (!is_water) continue;
            }
            // if !waters_only_flag, we process all residues

            // We need to iterate atoms and potentially delete some.
            // Deleting while iterating is dangerous, so we collect indices first.
            std::vector<size_t> indices_to_remove;

            if (remove_or_zero_occ_flag == TRIM_BY_MAP_DELETE) {
                for (size_t iat = 0; iat < residue.atoms.size(); iat++) {
                    const gemmi::Atom &at = residue.atoms[iat];
                    clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
                    if (coot::util::density_at_point(xmap, co) < map_level) {
                        indices_to_remove.push_back(iat);
                        n_changed++;
                    }
                }
            } else if (remove_or_zero_occ_flag == TRIM_BY_MAP_ZERO_OCC) {
                for (size_t iat = 0; iat < residue.atoms.size(); iat++) {
                    gemmi::Atom &at = residue.atoms[iat];
                    clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
                    if (coot::util::density_at_point(xmap, co) < map_level) {
                        at.occ = 0.0;
                        n_changed++;
                    }
                }
            }

            // Handle deletion: remove in reverse order to preserve indices
            if (!indices_to_remove.empty()) {
                for (auto it = indices_to_remove.rbegin(); it != indices_to_remove.rend(); ++it) {
                    residue.atoms.erase(residue.atoms.begin() + *it);
                }
            }
        }
    }

    return n_changed;
}

} // namespace util
} // namespace coot
