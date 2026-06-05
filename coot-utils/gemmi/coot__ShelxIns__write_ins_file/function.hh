#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/unitcell.hpp>
#include <clipper/clipper.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cctype>

namespace coot {

// reshelx: merge all chains into one chain "A", clear cell
inline gemmi::Structure reshelx_struct(const gemmi::Structure& st) {
  gemmi::Structure shelx_st;
  gemmi::Model model;
  gemmi::Chain chain;
  chain.name = "A";

  if (!st.models.empty()) {
    const gemmi::Model& src_model = st.models[0];
    for (const gemmi::Chain& src_chain : src_model.chains) {
      for (const gemmi::Residue& src_residue : src_chain.residues) {
        gemmi::Residue copy_residue = src_residue;
        copy_residue.atoms = src_residue.atoms;
        chain.residues.push_back(copy_residue);
      }
    }
  }

  model.chains.push_back(chain);
  shelx_st.models.push_back(model);

  shelx_st.cell.a = 0;
  shelx_st.cell.b = 0;
  shelx_st.cell.c = 0;
  shelx_st.cell.alpha = 0;
  shelx_st.cell.beta = 0;
  shelx_st.cell.gamma = 0;

  if (!st.spacegroup_hm.empty()) {
    shelx_st.spacegroup_hm = st.spacegroup_hm;
  }

  return shelx_st;
}

class ShelxIns {
public:
  clipper::Cell cell;
  int have_cell_flag = 0;

  bool try_assign_cell_gemmi(const gemmi::Structure& st) {
    if (!have_cell_flag) {
      double a = st.cell.a;
      double b = st.cell.b;
      double c = st.cell.c;
      double alpha_deg = st.cell.alpha;
      double beta_deg  = st.cell.beta;
      double gamma_deg = st.cell.gamma;

      clipper::Cell_descr cdr(a, b, c,
                              clipper::Util::d2rad(alpha_deg),
                              clipper::Util::d2rad(beta_deg),
                              clipper::Util::d2rad(gamma_deg));
      cell = clipper::Cell(cdr);
      have_cell_flag = 1;
    }
    return have_cell_flag;
  }

  bool mol_needs_shelx_transfer_gemmi(const gemmi::Structure& st) const {
    bool needs_unshelx = true;
    if (st.models.empty()) return true;
    const gemmi::Model& model = st.models[0];
    if (model.chains.size() > 1) {
      needs_unshelx = false;
    } else {
      for (const gemmi::Chain& chain : model.chains) {
        if (!chain.name.empty()) {
          needs_unshelx = false;
          break;
        }
      }
    }
    return needs_unshelx;
  }

  std::pair<int, std::string> write_ins_file_internal_gemmi(
      const gemmi::Model& model,
      const std::string& filename,
      bool mol_is_from_shelx_ins);

  std::pair<int, std::string> write_ins_file_gemmi(
      const gemmi::Structure& mol_in,
      const std::string& filename,
      bool mol_is_from_shelx_ins) {

    std::pair<int, std::string> r(-1, "");

    if (!have_cell_flag) {
      have_cell_flag = try_assign_cell_gemmi(mol_in);
    }

    if (mol_is_from_shelx_ins) {
      bool mol_needs_reshelx = mol_needs_shelx_transfer_gemmi(mol_in);
      if (mol_needs_reshelx) {
        gemmi::Structure mol = reshelx_struct(mol_in);
        if (mol.models.empty()) {
          return std::make_pair(-1, "ERROR:: no models after reshelx");
        }
        r = write_ins_file_internal_gemmi(mol.models[0], filename, true);
      } else {
        r = write_ins_file_internal_gemmi(mol_in.models[0], filename, true);
      }
    } else {
      r = write_ins_file_internal_gemmi(mol_in.models[0], filename, false);
    }

    return r;
  }
};

} // namespace coot