#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include <vector>
#include <map>
#include <string>

namespace coot {

/* ── gemmi-compatible overlap entry ── */
struct atom_overlap_t_gemmi {
    gemmi::Atom const* atom_1{nullptr};
    gemmi::Atom const* atom_2{nullptr};
    double overlap_volume{0.0};
};

/* ── gemmi-compatible result carrying overlaps data ── */
class atom_overlaps_container_t_gemmi_result {
public:
    std::vector<gemmi::Residue const*> neighbours;
    gemmi::Residue const* res_central{nullptr};
    gemmi::Model const* mol{nullptr};
    protein_geometry const* geom_p{nullptr};
    int imol_enc{0};
    double clash_spike_length{0.5};

    float score();
    void make_overlaps();
    double get_vdw_radius_atom(gemmi::Atom const* at);

private:
    std::vector<atom_overlap_t_gemmi> overlaps;
    std::map<std::string, double> type_to_vdw_radius_map;

    double get_overlap_volume(double dist, double r1, double r2) const;
};

/* ── factory — the ported constructor (free function in coot namespace) ── */
atom_overlaps_container_t_gemmi_result atom_overlaps_container_t_gemmi(
    gemmi::Residue const* res_central_in,
    gemmi::Residue const* neighbour,
    gemmi::Model const* mol_in,
    protein_geometry const* geom_p_in);

} // namespace coot