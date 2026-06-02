#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <cmath>
#include <tuple>

#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

// Our own enum to mirror coot::side_chain_densities::mode_t
enum side_chain_mode_t {
  MODE_GEN_USABLE_POINTS = 0,
  MODE_SAMPLE_FOR_DB     = 1,
  MODE_SAMPLE_FOR_RESIDUE = 2
};

// Our own density box type (avoids needing the original class / its private members)
struct density_box_t_gemmi {
   float *data;
   int    n_steps;
   int    size;
   double mean;
   double var;
   bool   around_ca;
   double ca_mean;
   double ca_var;
   double ca_sd;

   density_box_t_gemmi()
     : data(nullptr), n_steps(0), size(0), mean(0), var(0),
       around_ca(false), ca_mean(-1), ca_var(-1), ca_sd(-1) {}

   density_box_t_gemmi(float *d, int ns)
     : data(d), n_steps(ns), mean(0), var(0),
       around_ca(false), ca_mean(-1), ca_var(-1), ca_sd(-1) {
     if (d) {
       int n_per_side = 2 * ns + 1;
       size = n_per_side * n_per_side * n_per_side;
     } else {
       size = 0;
     }
   }

   bool empty() const { return (data == nullptr || size == 0); }

   void set_around_ca_stats(double m, double v, double s) {
     around_ca = true; ca_mean = m; ca_var = v; ca_sd = s;
   }
};

namespace coot {

/* ---- helpers (port of private side-chain-densities methods) ---- */

inline clipper::Coord_orth make_pt_in_grid_gemmi(
    int ix, int iy, int iz, float step_size,
    const std::vector<clipper::Coord_orth> &axes) {
   clipper::Coord_orth pt(0, 0, 0);
   int idx[3] = {ix, iy, iz};
   for (int i = 0; i < 3; ++i) {
      pt[0] += idx[i] * step_size * axes[i][0];
      pt[1] += idx[i] * step_size * axes[i][1];
      pt[2] += idx[i] * step_size * axes[i][2];
   }
   return pt;
}

inline bool is_close_to_atoms_gemmi(
    const std::vector<std::pair<double, clipper::Coord_orth>> &atom_positions,
    const clipper::Coord_orth &pt) {
   for (const auto &p : atom_positions) {
      clipper::Coord_orth d = pt - p.second;
      double dist = std::sqrt(d.lengthsq());
      if (dist < p.first) return true;
   }
   return false;
}

inline bool in_sphere_gemmi(
    const clipper::Coord_orth &pt,
    const clipper::Coord_orth &center,
    float radius) {
   clipper::Coord_orth d = pt - center;
   return std::sqrt(d.lengthsq()) <= radius;
}

// Trim whitespace so we can match CIF-style (unpadded) atom names
inline std::string trim_name(const std::string &s) {
   auto a = s.find_first_not_of(" \t\r\n");
   if (a == std::string::npos) return "";
   auto b = s.find_last_not_of(" \t\r\n");
   return s.substr(a, b - a + 1);
}

/* ---- the ported function ---- */

namespace side_chain_densities {

density_box_t_gemmi sample_map_gemmi(
    const gemmi::Residue *residue_this_p,
    const gemmi::Residue *residue_next_p,
    side_chain_mode_t mode,
    const clipper::Coord_orth &cb_pt,
    const std::vector<clipper::Coord_orth> &axes,
    const clipper::Xmap<float> &xmap,
    std::string gen_pts_file_name,
    float grid_box_radius,
    int n_steps,
    const std::set<int> &useable_grid_points);

} // namespace side_chain_densities
} // namespace coot