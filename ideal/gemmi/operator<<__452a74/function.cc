#include "function.hh"

namespace coot {

std::ostream &
operator_stream_insert_gemmi(std::ostream &s, const coot::phi_psi_t &pp) {
   s << pp.phi << " " << pp.psi;
   return s;
}

} // namespace coot