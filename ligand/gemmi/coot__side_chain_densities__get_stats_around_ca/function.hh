#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include "analysis/stats.hh"
#include "coot-utils/coot-map-utils.hh"
#include <tuple>
#include <vector>
#include <string>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

static inline clipper::Coord_orth
normalize_coord(const clipper::Coord_orth& v) {
    double len = std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    if (len == 0.0) return clipper::Coord_orth(0, 0, 0);
    return (static_cast<double>(1.0) / len) * v;
}

static inline std::vector<clipper::Coord_orth>
make_axes_gemmi(
    const clipper::Coord_orth& pt_ca,
    const clipper::Coord_orth& pt_cb,
    const clipper::Coord_orth& pt_c,
    const clipper::Coord_orth& pt_n
) {
    clipper::Coord_orth v_cb_to_ca = pt_ca - pt_cb;
    clipper::Coord_orth v_cb_ca = normalize_coord(v_cb_to_ca);
    clipper::Coord_orth v_n_to_c = pt_c - pt_n;
    clipper::Coord_orth v_n_c = normalize_coord(v_n_to_c);
    clipper::Coord_orth x_vec = v_cb_ca;
    clipper::Coord_orth z_vec = normalize_coord(clipper::Coord_orth(clipper::Coord_orth::cross(x_vec, v_n_c)));
    clipper::Coord_orth y_vec = normalize_coord(clipper::Coord_orth(clipper::Coord_orth::cross(z_vec, x_vec)));
    return {x_vec, y_vec, z_vec};
}

static inline clipper::Coord_orth
make_pt_in_grid_gemmi(int ix, int iy, int iz, float step_size,
                       const std::vector<clipper::Coord_orth>& axes) {
    return (static_cast<double>(ix) * static_cast<double>(step_size)) * axes[0]
         + (static_cast<double>(iy) * static_cast<double>(step_size)) * axes[1]
         + (static_cast<double>(iz) * static_cast<double>(step_size)) * axes[2];
}

std::tuple<double, double, double>
get_stats_around_ca_gemmi(
    const gemmi::Residue* residue,
    const std::vector<clipper::Coord_orth>& axes,
    float step_size,
    const clipper::Xmap<float>& xmap,
    int n_steps
) {
    double mean = 0;
    double var = 0;
    double mean_of_positives = 0;

    // Find CA atom
    const gemmi::Atom* ca_at = nullptr;
    for (const auto& at : residue->atoms) {
        std::string aname = at.name;
        size_t start = aname.find_first_not_of(" ");
        size_t end = aname.find_last_not_of(" ");
        std::string trimmed;
        if (start != std::string::npos) {
            trimmed = aname.substr(start, end - start + 1);
        }
        if (trimmed == "CA") {
            ca_at = &at;
            break;
        }
    }

    if (ca_at) {
        stats::single s;
        stats::single s_above_zero;
        clipper::Coord_orth ca_pos = coot::co_gemmi(ca_at);
        int ilim = n_steps * n_steps;

        for (int ix = -n_steps; ix <= n_steps; ix++) {
            for (int iy = -n_steps; iy <= n_steps; iy++) {
                for (int iz = -n_steps; iz <= n_steps; iz++) {
                    if ((ix * ix + iy * iy + iz * iz) <= ilim) {
                        clipper::Coord_orth pt_in_grid =
                            make_pt_in_grid_gemmi(ix, iy, iz, step_size, axes);
                        clipper::Coord_orth pt_grid_point = ca_pos + pt_in_grid;
                        float dv =
                            coot::util::density_at_point_by_linear_interpolation(xmap, pt_grid_point);
                        s.add(dv);
                        if (dv > 0.0)
                            s_above_zero.add(dv);
                    }
                }
            }
        }

        mean = s.mean();
        var = s.variance();
        mean_of_positives = s_above_zero.mean();

    } else {
        mean = -1;
        var = -1;
        mean_of_positives = -1;
    }

    return std::make_tuple(mean, var, mean_of_positives);
}

} // namespace coot