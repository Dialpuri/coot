#pragma once

#include <gemmi/model.hpp>

namespace coot {

class density_box_t {
public:
  density_box_t(float *density_box, gemmi::Residue *residue_p, int n_steps);
  density_box_t() { init(); }
  
  float *density_box;
  gemmi::Residue *residue_p;
  double mean;
  double mean_of_positives;
  double var;

  double mean_around_ca;
  double mean_of_positives_around_ca;
  double var_around_ca;

  bool is_weird;

  int n_steps;

  void init() {
    density_box = nullptr; residue_p = nullptr; n_steps = 0; mean=0; var = -1;
    mean_around_ca = 0; mean_of_positives_around_ca = 0;
    var_around_ca = -1;
    is_weird = false;
    mean_of_positives = 0;
  }
  
  bool empty() const { return (n_steps == 0); }
  
  int nnn() const {
    int n = 2 * n_steps + 1;
    return n * n * n;
  }
};

} // namespace coot