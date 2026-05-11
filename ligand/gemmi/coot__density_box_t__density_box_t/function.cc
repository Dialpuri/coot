#include "function.hh"

namespace coot {

density_box_t::density_box_t(float *density_box, gemmi::Residue *residue_p, int n_steps) {
  init();
  this->density_box = density_box;
  this->residue_p = residue_p;
  this->n_steps = n_steps;
}

} // namespace coot