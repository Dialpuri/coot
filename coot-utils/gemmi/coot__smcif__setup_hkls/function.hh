#pragma once

#include <string>
#include <vector>
#include <gemmi/cif.hpp>
#include "clipper/core/hkl_data.h"

namespace coot {
namespace smcif {

/// Gemmi-port of setup_hkls: reads a CIF file, extracts HKL indices
/// from _refln_index_h/k/l (or _pd_refln_index_h/k/l) loop.
/// Returns the list of HKLs found.
std::vector<clipper::HKL> setup_hkls_gemmi(const std::string& file_name);

}  // namespace smcif
}  // namespace coot