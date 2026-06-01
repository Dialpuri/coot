#include "function.hh"
#include <gemmi/read_cif.hpp>
#include <iostream>
#include <string>
#include <cmath>

namespace coot {
namespace smcif {

clipper::Resolution get_resolution_gemmi(const clipper::Cell &cell,
                                         const std::string &file_name) {
  clipper::ftype slim = 0.0;

  gemmi::cif::Document doc;
  try {
    doc = gemmi::read_cif_gz(file_name);
  } catch (const std::exception&) {
    std::cout << "WARNING:: Error reading small-molecule cif \"" << file_name << "\"" << std::endl;
    double reso_A = 1.0 / std::sqrt(slim);
    return clipper::Resolution(reso_A);
  }

  auto& blocks = doc.blocks;
  if (blocks.empty()) {
    double reso_A = 1.0 / std::sqrt(slim);
    return clipper::Resolution(reso_A);
  }

  auto& block = blocks[0];

  // Try _refln_index_* tags first, then _pd_refln_index_*
  gemmi::cif::Table tab_refln = block.find({"_refln_index_h", "_refln_index_k", "_refln_index_l"});
  gemmi::cif::Table tab_pd = block.find({"_pd_refln_index_h", "_pd_refln_index_k", "_pd_refln_index_l"});

  // Pick whichever has data
  gemmi::cif::Table* tab = nullptr;
  if (tab_refln.size() > 0) {
    tab = &tab_refln;
  } else if (tab_pd.size() > 0) {
    tab = &tab_pd;
  }

  if (tab) {
    for (auto row : *tab) {
      std::string h_str = row.str(0);
      if (h_str.empty()) continue;
      std::string k_str = row.str(1);
      if (k_str.empty()) continue;
      std::string l_str = row.str(2);
      if (l_str.empty()) continue;

      int h = std::stoi(h_str);
      int k = std::stoi(k_str);
      int l = std::stoi(l_str);
      clipper::HKL hkl(h, k, l);
      double reso = hkl.invresolsq(cell);
      slim = clipper::Util::max(slim, reso);
    }
  }

  double reso_A = 1.0 / std::sqrt(slim);
  return clipper::Resolution(reso_A);
}

} // namespace smcif
} // namespace coot