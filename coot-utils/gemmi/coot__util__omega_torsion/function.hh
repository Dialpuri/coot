#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot { namespace util {

std::pair<bool, double>
omega_torsion_gemmi(const gemmi::Residue* C_residue,
                    const gemmi::Residue* N_residue,
                    const std::string& altconf);

}} // namespace coot::util