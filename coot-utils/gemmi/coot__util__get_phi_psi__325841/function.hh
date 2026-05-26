#pragma once
#include <utility>
#include <string>
#include <gemmi/model.hpp>

namespace coot { namespace util {

struct phi_psi_with_residues_t {
    double phi_;
    double psi_;
    const gemmi::Residue* residue_prev;
    const gemmi::Residue* residue_this;
    const gemmi::Residue* residue_next;

    phi_psi_with_residues_t()
        : phi_(0.0), psi_(0.0),
          residue_prev(nullptr), residue_this(nullptr), residue_next(nullptr) {}
};

// Port of coot::util::get_phi_psi(mmdb::PResidue*) to gemmi
// Takes array of 3 consecutive gemmi residues and returns phi/psi dihedrals
std::pair<bool, phi_psi_with_residues_t>
get_phi_psi_gemmi(const gemmi::Residue* residue_0,
                  const gemmi::Residue* residue_1,
                  const gemmi::Residue* residue_2);

}} // namespace coot::util