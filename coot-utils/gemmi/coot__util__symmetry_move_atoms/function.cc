#include "function.hh"
#include "coot-utils/coot-coord-utils.hh"

namespace coot {
namespace util {

std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>>
symmetry_move_atoms_gemmi(
    const std::vector<clipper::Coord_orth>& protein_coords,
    const std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>>& water_atoms_in,
    clipper::Cell cell,
    clipper::Spacegroup spacegroup) {

  std::vector<double> current_best_dist(water_atoms_in.size(), -1);

  clipper::Coord_frac pre_shift_frac =
      shift_to_origin(protein_coords, cell, spacegroup);
  clipper::Coord_orth pre_shift_orth = pre_shift_frac.coord_orth(cell);

  // Let's not move atoms if the best distance is not better than the current
  // best distance so let's make a note of those now.
  std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>> water_atoms =
      water_atoms_in;
  for (unsigned int iat = 0; iat < water_atoms_in.size(); iat++) {
    for (unsigned int iprot = 0; iprot < protein_coords.size(); iprot++) {
      double d = clipper::Coord_orth::length(
          water_atoms_in[iat].second, protein_coords[iprot]);
      if (current_best_dist[iat] < 0) {
        current_best_dist[iat] = d;
      } else {
        if (d < current_best_dist[iat])
          current_best_dist[iat] = d;
      }
    }
  }

  // create shifted protein coords
  std::vector<clipper::Coord_orth> protein_coords_origin_shifted(
      protein_coords.size());
  for (unsigned int ip = 0; ip < protein_coords.size(); ip++) {
    protein_coords_origin_shifted[ip] = protein_coords[ip] + pre_shift_orth;
  }

  // Do the cell shift search
  int n = spacegroup.num_symops();
  clipper::Coord_frac cell_shift;
  for (unsigned int iw = 0; iw < water_atoms.size(); iw++) {
    clipper::Coord_orth water_pos_pre(water_atoms[iw].second);
    clipper::Coord_orth water_pos =
        translate_close_to_origin(water_pos_pre, cell);

    double d_best = 99999999.9;
    // The compiler doesn't like rtop_best being used below
    // without being initialized properly here.
    clipper::RTop_orth rtop_best(
        clipper::Mat33<double>(1, 0, 0, 0, 1, 0, 0, 0, 1),
        clipper::Coord_orth(0, 0, 0));
    bool improved = false;

    for (int isym = 0; isym < n; isym++) {
      for (int x_shift = -1; x_shift < 2; x_shift++) {
        for (int y_shift = -1; y_shift < 2; y_shift++) {
          for (int z_shift = -1; z_shift < 2; z_shift++) {
            cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
            clipper::RTop_orth orthop =
                clipper::RTop_frac(
                    spacegroup.symop(isym).rot(),
                    spacegroup.symop(isym).trn() + cell_shift)
                    .rtop_orth(cell);
            clipper::Coord_orth t_point = water_pos.transform(orthop);
            double t_dist =
                min_dist_to_points(t_point, protein_coords_origin_shifted);
            if (t_dist < d_best) {
              if (t_dist < current_best_dist[iw]) {
                d_best = t_dist;
                rtop_best = orthop;
                improved = true;
              }
            }
          }
        }
      }
    }

    if (improved) {
      // Apply the transformation then.
      clipper::Coord_orth t_point = water_pos.transform(rtop_best);
      water_atoms[iw].second = t_point - pre_shift_orth;
    } else {
      water_atoms[iw].first = nullptr; // NULL, don't move it.
    }
  }
  return water_atoms;
}

} // namespace util
} // namespace coot