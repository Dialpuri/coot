#include "function.hh"
#include "geometry/protein-geometry.hh"
#include <cmath>

namespace coot::hole {

double sphere_size_gemmi(
    const gemmi::Model& model,
    const std::vector<std::string>& chain_names,
    const gemmi::Vec3& pt,
    const coot::protein_geometry& geom
) {
    double largest_possible_sphere = 99999.0;
    bool was_set = false;

    for (const auto& chain : model.chains) {
        // If chain_names is non-empty, only process matching chains
        if (!chain_names.empty()) {
            bool found = false;
            for (const auto& cn : chain_names) {
                if (chain.name == cn) {
                    found = true;
                    break;
                }
            }
            if (!found) continue;
        }

        for (const auto& residue : chain.residues) {
            std::string res_name = residue.name;

            for (const auto& atom : residue.atoms) {
                // Compute distance from pt to atom position
                double dx = atom.pos.x - pt.x;
                double dy = atom.pos.y - pt.y;
                double dz = atom.pos.z - pt.z;
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                // Get VDW radius from protein_geometry
                double vdw_radius = geom.get_vdw_radius(
                    atom.name, res_name, 0, false);

                // Fallback: if get_vdw_radius returns <= 0, use element-based default
                if (vdw_radius <= 0) {
                    vdw_radius = 1.7;
                    std::string ele = atom.element.name();
                    if (ele == "N")  vdw_radius = 1.55;
                    if (ele == "O")  vdw_radius = 1.52;
                    if (ele == "H")  vdw_radius = 1.2;
                }

                double r = dist - vdw_radius;
                if (r < largest_possible_sphere) {
                    largest_possible_sphere = r;
                    if (!was_set)
                        was_set = true;
                }
            }
        }
    }

    if (was_set)
        return largest_possible_sphere;
    else
        return -1.0; // no atoms selected
}

} // namespace coot::hole