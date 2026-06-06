#pragma once

#include <string>
#include <gemmi/model.hpp>

// Co-dependencies
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__primitive_chi_angles__get_chi_angles/gemmi/function.hh"

namespace rama_plot {
class phi_psi_t {
public:
    double phi; // in degrees
    double psi; // in degrees
    phi_psi_t(double phi_val = 0, double psi_val = 0) : phi(phi_val), psi(psi_val) {}
};
}

namespace coot {
class molecule_t {
public:
    static std::string get_torsions_for_residues_in_chain_as_json_gemmi(
        gemmi::Structure& structure, const std::string& chain_id);
};
} // namespace coot