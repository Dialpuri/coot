#include "function.hh"
#include <cmath>

namespace coot {

/* ── Factory (ported constructor) ── */
atom_overlaps_container_t_gemmi_result atom_overlaps_container_t_gemmi(
    gemmi::Residue const* res_central_in,
    gemmi::Residue const* neighbour,
    gemmi::Model const* mol_in,
    protein_geometry const* geom_p_in)
{
    atom_overlaps_container_t_gemmi_result aoc;
    aoc.geom_p            = geom_p_in;
    aoc.res_central       = res_central_in;
    aoc.mol               = mol_in;
    aoc.imol_enc          = protein_geometry::IMOL_ENC_ANY;
    aoc.clash_spike_length = 0.5;
    aoc.neighbours.push_back(neighbour);
    aoc.make_overlaps();
    return aoc;
}

/* ── Overlap helpers ── */
double atom_overlaps_container_t_gemmi_result::get_vdw_radius_atom(gemmi::Atom const* at)
{
    auto it = type_to_vdw_radius_map.find(at->element.name());
    if (it != type_to_vdw_radius_map.end()) {
        return it->second;
    }

    double vdwr = 2.0; // default
    if (geom_p && res_central) {
        vdwr = geom_p->get_vdw_radius(at->element.name(),
                                       res_central->name,
                                       imol_enc,
                                       true);
    }

    type_to_vdw_radius_map[at->element.name()] = vdwr;
    return vdwr;
}

double atom_overlaps_container_t_gemmi_result::get_overlap_volume(
    double dist, double r1, double r2) const
{
    double sum = r1 + r2;
    double dif = r1 - r2;
    if (dist >= sum || dist <= dif) {
        return 0.0;
    }
    double v = (M_PI / 12.0) * (sum + dist) * (sum - dist) * (sum + dist);
    return std::max(v, 0.0);
}

void atom_overlaps_container_t_gemmi_result::make_overlaps()
{
    overlaps.clear();
    type_to_vdw_radius_map.clear();

    if (!res_central) return;
    if (neighbours.empty()) return;

    for (auto* nbr : neighbours) {
        if (!nbr) continue;
        for (auto& atom_c : res_central->atoms) {
            for (auto& atom_n : nbr->atoms) {
                if (&atom_c == &atom_n) continue;
                double dx = atom_c.pos.x - atom_n.pos.x;
                double dy = atom_c.pos.y - atom_n.pos.y;
                double dz = atom_c.pos.z - atom_n.pos.z;
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < clash_spike_length) {
                    double r1 = get_vdw_radius_atom(&atom_c);
                    double r2 = get_vdw_radius_atom(&atom_n);
                    double ov = get_overlap_volume(dist, r1, r2);
                    if (ov > 0.0) {
                        atom_overlap_t_gemmi o;
                        o.atom_1         = &atom_c;
                        o.atom_2         = &atom_n;
                        o.overlap_volume = ov;
                        overlaps.push_back(o);
                    }
                }
            }
        }
    }
}

float atom_overlaps_container_t_gemmi_result::score()
{
    float s = 0.0f;
    unsigned int nos = static_cast<unsigned int>(overlaps.size());
    if (nos > 0) {
        for (unsigned int i = 0; i < nos; i++) {
            s += static_cast<float>(overlaps[i].overlap_volume);
        }
        s /= static_cast<float>(nos);
        s *= 1000.0f;
    }
    return s;
}

} // namespace coot