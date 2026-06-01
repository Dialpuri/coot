#pragma once

#include <sstream>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>

namespace coot {

// Port of coot::molecule_t::molecule_to_mmCIF_string
// Uses gemmi::make_mmcif_document + gemmi::cif::write_cif_to_stream to convert a gemmi::Structure
// to an mmCIF string.
inline std::string molecule_to_mmCIF_string_gemmi(const gemmi::Structure& st) {
  if (st.models.empty()) {
    return std::string();
  }
  auto doc = gemmi::make_mmcif_document(st);
  std::ostringstream oss;
  gemmi::cif::write_cif_to_stream(oss, doc);
  return oss.str();
}

} // namespace coot