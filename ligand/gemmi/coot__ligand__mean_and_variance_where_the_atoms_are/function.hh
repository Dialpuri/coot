#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include <clipper/core/map_interp.h>
#include "coot/analysis/stats.hh"
#include "utils/coot-utils.hh"
#include <utility>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>

namespace coot {
namespace ligand {

// Replicates density_at_point from the original coot::ligand class
inline float density_at_point_gemmi(const clipper::Coord_orth &atom_pos, const clipper::Xmap<float> &search_map) {
   clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(search_map.cell());
   float dv = search_map.interp<clipper::Interp_cubic>(atom_pos_frc);
   return dv;
}

inline std::pair<float, float> mean_and_variance_where_the_atoms_are_gemmi(
    const gemmi::Model &model,
    const clipper::Xmap<float> &xmap)
{
  std::pair<float, float> r(0, 0);
  unsigned int n_test_points = 100;
  std::vector<clipper::Coord_orth> test_points;

  unsigned int n_molecule_atoms = 0;

  // First pass: count non-H, non-terminator atoms (skip HOH)
  for (const gemmi::Chain &chain : model.chains) {
    for (const gemmi::Residue &residue : chain.residues) {
      std::string rn = residue.name;
      if (rn != "HOH") {
        for (const gemmi::Atom &atom : residue.atoms) {
          // Check for terminators (residue names ending in "XT")
          bool is_term = (rn.size() >= 2 && rn.substr(rn.size() - 2) == "XT");
          if (!is_term) {
            if (!atom.element.is_hydrogen()) {
              n_molecule_atoms++;
            }
          }
        }
      }
    }
  }

  // Second pass: collect test points
  if (n_molecule_atoms > n_test_points) {
    float rmi = 1.0f / static_cast<float>(RAND_MAX);
    float crit_val = static_cast<float>(n_test_points) / static_cast<float>(n_molecule_atoms);
    for (const gemmi::Chain &chain : model.chains) {
      for (const gemmi::Residue &residue : chain.residues) {
        std::string rn = residue.name;
        if (rn != "HOH") {
          for (const gemmi::Atom &atom : residue.atoms) {
            bool is_term = (rn.size() >= 2 && rn.substr(rn.size() - 2) == "XT");
            if (!is_term) {
              if (!atom.element.is_hydrogen()) {
                float f = coot::util::random() * rmi;
                if (f < crit_val) {
                  clipper::Coord_orth c(atom.pos.x, atom.pos.y, atom.pos.z);
                  test_points.push_back(c);
                }
              }
            }
          }
        }
      }
    }
  } else {
    for (const gemmi::Chain &chain : model.chains) {
      for (const gemmi::Residue &residue : chain.residues) {
        std::string rn = residue.name;
        if (rn != "HOH") {
          for (const gemmi::Atom &atom : residue.atoms) {
            bool is_term = (rn.size() >= 2 && rn.substr(rn.size() - 2) == "XT");
            if (!is_term) {
              if (!atom.element.is_hydrogen()) {
                clipper::Coord_orth c(atom.pos.x, atom.pos.y, atom.pos.z);
                test_points.push_back(c);
              }
            }
          }
        }
      }
    }
  }

  if (!test_points.empty()) {
    coot::stats::single s;
    for (std::size_t i = 0; i < test_points.size(); i++) {
      float val = density_at_point_gemmi(test_points[i], xmap);
      s.add(val);
    }
    float m = s.mean();
    float sd = std::sqrt(s.variance());
    return std::pair<float, float>(m, sd);
  }

  return r;
}

} // namespace ligand
} // namespace coot