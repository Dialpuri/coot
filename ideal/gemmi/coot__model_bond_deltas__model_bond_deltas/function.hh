#pragma once

#include <map>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

class xyz_deltas_t {
public:
  std::vector<double> x[3];
  std::vector<double> deltas; // not resolved along axes
  unsigned int n;
  xyz_deltas_t() { n = 0; }
  void add(double delta_in, const clipper::Coord_orth &d) {
    deltas.push_back(delta_in);
    x[0].push_back(d.x());
    x[1].push_back(d.y());
    x[2].push_back(d.z());
    ++n;
  }
  unsigned int size() const { return n; }
};

class model_bond_deltas_gemmi
    : public std::map<std::string, std::map<int, std::map<int, xyz_deltas_t>>> {
public:
  gemmi::Model *mol;
  int imol;
  coot::protein_geometry *geom_p;

  model_bond_deltas_gemmi()
    : mol(nullptr),
      imol(coot::protein_geometry::IMOL_ENC_UNSET),
      geom_p(nullptr)
  {
  }

  model_bond_deltas_gemmi(gemmi::Model *mol_in,
                          int imol_in,
                          coot::protein_geometry *geom_p_in)
    : mol(mol_in),
      imol(imol_in),
      geom_p(geom_p_in)
  {
  }

  xyz_deltas_t xyzd;
  void add(double delta_in, const clipper::Coord_orth &d) {
    xyzd.add(delta_in, d);
  }
  unsigned int size() const { return xyzd.size(); }
};

} // namespace coot