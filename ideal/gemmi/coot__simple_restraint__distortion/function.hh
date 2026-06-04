#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <utility>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__angle/gemmi/function.hh"

// Enum is at global scope (matches coot/ideal/simple-restraint.hh)
enum restraint_type_t {
  BOND_RESTRAINT = 1,
  ANGLE_RESTRAINT,
  CHIRAL_VOLUME_RESTRAINT,
  NON_BONDED_CONTACT_RESTRAINT,
  GEMAN_MCCLURE_DISTANCE_RESTRAINT
};

namespace coot {

class simple_restraint {
public:
  restraint_type_t restraint_type = BOND_RESTRAINT;
  int atom_index_1    = -1;
  int atom_index_2    = -1;
  int atom_index_3    = -1;
  int atom_index_centre = -1;
  double target_value        = 0.0;
  double target_chiral_volume = 0.0;
  double sigma = 1.0;

  std::pair<double, double>
  distortion_gemmi(std::vector<gemmi::Atom*> atoms, const double &lj_epsilon) const {

    std::pair<double, double> distortion_pair(-1.0, -1.0);

    if (restraint_type == CHIRAL_VOLUME_RESTRAINT) {
      gemmi::Atom *at_c = atoms[atom_index_centre];
      gemmi::Atom *at_1 = atoms[atom_index_1];
      gemmi::Atom *at_2 = atoms[atom_index_2];
      gemmi::Atom *at_3 = atoms[atom_index_3];

      clipper::Coord_orth centre = coot::co_gemmi(at_c);
      clipper::Coord_orth a1 = coot::co_gemmi(at_1);
      clipper::Coord_orth a2 = coot::co_gemmi(at_2);
      clipper::Coord_orth a3 = coot::co_gemmi(at_3);

      clipper::Coord_orth a = a1 - centre;
      clipper::Coord_orth b = a2 - centre;
      clipper::Coord_orth c = a3 - centre;

      double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b, c));
      double delta = cv - target_chiral_volume;
      distortion_pair.first  = delta * delta / (sigma * sigma);
      distortion_pair.second = delta;
    }

    if (restraint_type == BOND_RESTRAINT) {
      gemmi::Atom *at_1 = atoms[atom_index_1];
      gemmi::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
        clipper::Coord_orth p1 = coot::co_gemmi(at_1);
        clipper::Coord_orth p2 = coot::co_gemmi(at_2);
        double d = std::sqrt((p2 - p1).lengthsq());
        double delta = d - target_value;
        double z = delta / sigma;
        double pen_score = z * z;
        distortion_pair.first  = pen_score;
        distortion_pair.second = delta;
      }
    }

    if (restraint_type == GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
      gemmi::Atom *at_1 = atoms[atom_index_1];
      gemmi::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
        clipper::Coord_orth p1 = coot::co_gemmi(at_1);
        clipper::Coord_orth p2 = coot::co_gemmi(at_2);
        double d = std::sqrt((p2 - p1).lengthsq());
        double alpha = 0.01;
        double delta = d - target_value;
        double z = delta / sigma;
        double distortion = z * z / (1 + alpha * z * z);
        distortion_pair.first  = distortion;
        distortion_pair.second = delta;
      }
    }

    if (restraint_type == ANGLE_RESTRAINT) {
      gemmi::Atom *at_1 = atoms[atom_index_1];
      gemmi::Atom *at_2 = atoms[atom_index_2];
      gemmi::Atom *at_3 = atoms[atom_index_3];
      if (at_1 && at_2 && at_3) {
        double angle_deg = coot::angle_gemmi(at_1, at_2, at_3);
        double delta = angle_deg - target_value;
        double z = delta / sigma;
        double distortion = z * z;
        distortion_pair.first  = distortion;
        distortion_pair.second = delta;
      }
    }

    if (restraint_type == NON_BONDED_CONTACT_RESTRAINT) {
      gemmi::Atom *at_1 = atoms[atom_index_1];
      gemmi::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
        clipper::Coord_orth p1 = coot::co_gemmi(at_1);
        clipper::Coord_orth p2 = coot::co_gemmi(at_2);

        double dist_sq = clipper::Coord_orth(p2 - p1).lengthsq();
        double dist = std::sqrt(dist_sq);
        double dist_delta = dist - target_value;

        float V_lj = 0;

        double lj_sigma = target_value;
        double lj_r_min = 1.122462048309373 * lj_sigma;

        if (dist_sq < 0.81) dist_sq = 0.81;
        double alpha_sqrd  = lj_r_min * lj_r_min / dist_sq;
        double alpha_up_6  = alpha_sqrd * alpha_sqrd * alpha_sqrd;
        double alpha_up_12 = alpha_up_6 * alpha_up_6;
        V_lj = static_cast<float>(lj_epsilon * (alpha_up_12 - 2.0 * alpha_up_6));

        double Vlj_at_rmax = -0.016316891136 * lj_epsilon;
        V_lj += static_cast<float>(Vlj_at_rmax);

        distortion_pair.first  = static_cast<double>(V_lj);
        distortion_pair.second = dist_delta;
      }
    }

    return distortion_pair;
  }
};

} // namespace coot