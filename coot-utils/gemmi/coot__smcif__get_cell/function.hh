#pragma once

#include <gemmi/read_cif.hpp>
#include "clipper/core/cell.h"
#include "utils/coot-utils.hh"
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <sstream>

namespace coot {
namespace smcif {

clipper::Cell get_cell_gemmi(const std::string& filename) {
  std::string cell_a, cell_b, cell_c;
  std::string cell_alpha, cell_beta, cell_gamma;

  std::ifstream ifs(filename);
  std::string line;
  while (std::getline(ifs, line)) {
    // Skip empty lines and lines starting with '#'
    if (line.empty() || line[0] == '#') continue;

    if (line == "data_") {
      // Check next lines for the data block name
      continue;
    }
    if (line.rfind("data_", 0) == 0) continue;

    // Parse _tag value lines
    if (line[0] == '_') {
      size_t pos = line.find_first_of(" \t");
      if (pos != std::string::npos) {
        std::string tag = line.substr(0, pos);
        std::string val = line.substr(pos + 1);
        // Trim leading whitespace from value
        size_t vstart = val.find_first_not_of(" \t");
        if (vstart != std::string::npos)
          val = val.substr(vstart);
        // Trim trailing whitespace
        size_t vend = val.find_last_not_of(" \t\r\n");
        if (vend != std::string::npos)
          val = val.substr(0, vend + 1);

        if (tag == "_cell_length_a") cell_a = val;
        else if (tag == "_cell_length_b") cell_b = val;
        else if (tag == "_cell_length_c") cell_c = val;
        else if (tag == "_cell_angle_alpha") cell_alpha = val;
        else if (tag == "_cell_angle_beta") cell_beta = val;
        else if (tag == "_cell_angle_gamma") cell_gamma = val;
      }
    }
  }

  clipper::Cell cell;

  if (!cell_a.empty() && !cell_b.empty() && !cell_c.empty() &&
      !cell_alpha.empty() && !cell_beta.empty() && !cell_gamma.empty()) {
    std::vector<std::string> a_v     = coot::util::split_string_no_blanks(cell_a, "(");
    std::vector<std::string> b_v     = coot::util::split_string_no_blanks(cell_b, "(");
    std::vector<std::string> c_v     = coot::util::split_string_no_blanks(cell_c, "(");
    std::vector<std::string> alpha_v = coot::util::split_string_no_blanks(cell_alpha, "(");
    std::vector<std::string> beta_v  = coot::util::split_string_no_blanks(cell_beta,  "(");
    std::vector<std::string> gamma_v = coot::util::split_string_no_blanks(cell_gamma, "(");

    double a     = coot::util::string_to_float(a_v[0]);
    double b     = coot::util::string_to_float(b_v[0]);
    double c     = coot::util::string_to_float(c_v[0]);
    double alpha = coot::util::string_to_float(alpha_v[0]);
    double beta  = coot::util::string_to_float(beta_v[0]);
    double gamma = coot::util::string_to_float(gamma_v[0]);

    clipper::Cell_descr cell_descr(a, b, c,
                                   clipper::Util::d2rad(alpha),
                                   clipper::Util::d2rad(beta),
                                   clipper::Util::d2rad(gamma));
    cell.init(cell_descr);
  } else {
    throw std::runtime_error("smcif: failed to get cell");
  }

  return cell;
}

} // namespace smcif
} // namespace coot