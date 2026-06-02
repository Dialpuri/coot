#pragma once

#include <ostream>

namespace coot {

struct phi_psi_t {
    float phi;
    float psi;
    float tau;

    phi_psi_t(float a, float b) : phi(a), psi(b), tau(0.0f) {}
    phi_psi_t(float a, float b, float c) : phi(a), psi(b), tau(c) {}
};

std::ostream &
operator_stream_insert_gemmi(std::ostream &s, const coot::phi_psi_t &pp);

} // namespace coot