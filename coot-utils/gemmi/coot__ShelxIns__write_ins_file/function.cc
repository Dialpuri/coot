#include "function.hh"
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

namespace coot {

static std::string make_atom_name(const gemmi::Atom& atom) {
  return std::string(1, atom.altloc) + atom.name;
}

static int round_to_4_dp(double d) {
  return static_cast<int>(std::round(d * 10000.0));
}

static std::string format_double_5_4(double v) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(4) << std::setw(9) << std::right << v;
  return oss.str();
}

static bool cell_is_empty(const clipper::Cell& cell) {
  clipper::Cell_descr cd = cell.descr();
  return (cd.a() == 0 && cd.b() == 0 && cd.c() == 0);
}

static std::string synthetic_atom_line(const gemmi::Atom& atom, const clipper::Cell& cell) {
  double x = round_to_4_dp(atom.pos.x) / 10000.0;
  double y = round_to_4_dp(atom.pos.y) / 10000.0;
  double z = round_to_4_dp(atom.pos.z) / 10000.0;

  std::string atom_name = make_atom_name(atom);

  if (cell_is_empty(cell)) {
    std::ostringstream oss;
    oss << " " << atom_name << " " << format_double_5_4(x) << format_double_5_4(y)
        << format_double_5_4(z) << std::setw(6) << std::fixed << std::setprecision(2) << atom.occ
        << std::setw(7) << std::fixed << std::setprecision(2) << atom.b_iso;
    return oss.str();
  } else {
    clipper::Cell_descr cd = cell.descr();
    gemmi::UnitCell gc(cd.a(), cd.b(), cd.c(),
                       cd.alpha_deg(), cd.beta_deg(), cd.gamma_deg());
    gemmi::Fractional frac = gc.fractionalize(gemmi::Position(x, y, z));
    std::ostringstream oss;
    oss << "F " << atom_name << " " << format_double_5_4(frac.x) << format_double_5_4(frac.y)
        << format_double_5_4(frac.z) << std::setw(6) << std::fixed << std::setprecision(2) << atom.occ
        << std::setw(7) << std::fixed << std::setprecision(2) << atom.b_iso;
    return oss.str();
  }
}

static std::string residue_numbering(const gemmi::Residue& res, bool is_het) {
  std::ostringstream oss;
  int seq = res.seqid.num.value;
  char icode = res.seqid.icode;
  if (is_het) {
    oss << " " << std::setw(6) << std::right << seq;
  } else {
    if (icode != ' ') {
      oss << " " << std::setw(2) << std::right << seq << " " << icode;
    } else {
      oss << std::setw(3) << std::right << seq << "   ";
    }
  }
  return oss.str();
}

static std::string write_synthetic_hydrogen_bonds(
    const gemmi::Model& model,
    std::map<std::string, std::vector<gemmi::Residue*>>& res_per_chain) {
  std::string af_lines;

  auto find_res = [&res_per_chain](const gemmi::Atom& atom, int offset) -> gemmi::Residue* {
    for (auto& kv : res_per_chain) {
      for (auto* rp : kv.second) {
        for (const auto& a : rp->atoms) {
          if (&a == &atom) {
            if (offset == 0) return rp;
            auto& vec = kv.second;
            int idx = -1;
            for (int i = 0; i < static_cast<int>(vec.size()); i++) {
              if (vec[i] == rp) { idx = i; break; }
            }
            int new_idx = idx + offset;
            if (new_idx >= 0 && new_idx < static_cast<int>(vec.size())) return vec[new_idx];
            return nullptr;
          }
        }
      }
    }
    return nullptr;
  };

  for (const auto& chain : model.chains) {
    auto rit = res_per_chain.find(chain.name);
    if (rit == res_per_chain.end()) continue;
    const auto& residues = rit->second;

    for (int i = 0; i < static_cast<int>(residues.size()); i++) {
      const auto& r = residues[i];
      if (r->seqid.num.value > 0) {
        for (const auto& atom : r->atoms) {
          if (atom.element.name() == std::string("O")) {
            std::string donor;
            gemmi::Residue* donor_res = nullptr;
            if (r->name == "TYR") { donor = "OH"; donor_res = find_res(atom, 0); }
            else if (r->name == "PHE") { donor = "CZ"; donor_res = find_res(atom, 0); }
            else if (r->name == "ASP" && atom.name == "OD1") { donor = "OD2"; donor_res = find_res(atom, 0); }
            else if (r->name == "GLU" && atom.name == "OE1") { donor = "OE2"; donor_res = find_res(atom, 0); }
            else if (r->name == "GLN" && atom.name == "OE1") { donor = "NE2"; donor_res = find_res(atom, 0); }
            else if (r->name == "ASN" && atom.name == "OD1") { donor = "ND2"; donor_res = find_res(atom, 0); }
            else if (r->name == "GLU" && atom.name == "OE2") { donor = "OE1"; donor_res = find_res(atom, 0); }
            else if (r->name == "ASP" && atom.name == "OD2") { donor = "OD1"; donor_res = find_res(atom, 0); }
            else if (r->name == "HIS" && atom.name == "NE2") { donor = "ND1"; donor_res = find_res(atom, 0); }
            else if (r->name == "HIS" && atom.name == "ND1") { donor = "NE2"; donor_res = find_res(atom, 0); }

            if (!donor.empty() && donor_res) {
              for (const auto& a2 : donor_res->atoms) {
                if (a2.name == donor) {
                  std::ostringstream oss;
                  oss << "AFIX 147 " << make_atom_name(a2) << " " << chain.name << " "
                      << std::setw(2) << std::right << r->seqid.num.value << std::endl;
                  af_lines += oss.str();
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  return af_lines;
}

static void write_synthetic_pre_atom_lines_gemmi(
    const gemmi::Model& model,
    const clipper::Cell& cell,
    const std::vector<std::string>& sfac,
    std::ofstream& f) {

  f << "TITL PDB->ins" << std::endl;
  f << "CELL 0" << std::endl;
  f << "LATT 1" << std::endl;

  clipper::Cell_descr cd = cell.descr();
  if (cell_is_empty(cell)) {
    f << "SFAC " << sfac[0];
    for (size_t i = 1; i < sfac.size(); i++) f << " " << sfac[i];
    f << std::endl;
    f << std::endl;
  } else {
    f << std::fixed << std::setprecision(2) << std::setw(12) << cd.a()
      << std::setw(12) << cd.b() << std::setw(12) << cd.c()
      << std::setw(12) << cd.alpha_deg()
      << std::setw(12) << cd.beta_deg()
      << std::setw(12) << cd.gamma_deg() << std::endl;

    f << "ZERR ";
    f << std::setw(12) << 10 << std::setw(12) << 10 << std::setw(12) << 10
      << std::setw(12) << 10 << std::setw(12) << 10 << std::setw(12) << 10 << std::endl;

    f << "ZFIX ";
    f << std::setw(12) << 100 << std::setw(12) << 100 << std::setw(12) << 100
      << std::setw(12) << 100 << std::setw(12) << 100 << std::setw(12) << 100 << std::endl;

    f << "SFAC " << sfac[0];
    for (size_t i = 1; i < sfac.size(); i++) f << " " << sfac[i];
    f << std::endl;

    f << std::endl;
  }
}

std::pair<int, std::string> ShelxIns::write_ins_file_internal_gemmi(
    const gemmi::Model& model,
    const std::string& filename,
    bool mol_is_from_shelx_ins) {

  int n_residues = 0;
  std::ofstream f(filename.c_str());

  std::vector<std::string> sfac;
  std::set<std::string> sfac_set;
  std::map<std::string, std::vector<gemmi::Residue*>> res_per_chain;
  std::vector<gemmi::Residue*> res_with_hetflag;

  for (const gemmi::Chain& chain : model.chains) {
    for (const auto& residue : chain.residues) {
      n_residues++;
      res_per_chain[chain.name].push_back(const_cast<gemmi::Residue*>(&residue));
      for (const auto& atom : residue.atoms) {
        sfac_set.insert(atom.element.name());
      }
      if (residue.het_flag != ' ') {
        res_with_hetflag.push_back(const_cast<gemmi::Residue*>(&residue));
      }
    }
  }

  for (const auto& s : sfac_set) sfac.push_back(s);

  write_synthetic_pre_atom_lines_gemmi(model, cell, sfac, f);

  if (!mol_is_from_shelx_ins) {
    f << "HTAB ?" << std::endl;
    std::string af_lines = write_synthetic_hydrogen_bonds(model, res_per_chain);
    f << af_lines;
  }

  for (const gemmi::Chain& chain : model.chains) {
    for (const auto& residue : chain.residues) {
      bool is_het = (residue.het_flag != ' ');
      if (!is_het) {
        f << "C " << chain.name << " " << residue.name << " "
          << residue_numbering(residue, false) << std::endl;
      }
      for (const auto& atom : residue.atoms) {
        f << synthetic_atom_line(atom, cell) << std::endl;
      }
      if (is_het && res_with_hetflag.empty()) {
        f << "END" << std::endl;
        f << std::endl;
      }
    }
  }

  if (res_with_hetflag.empty()) {
    f << "END" << std::endl;
    f << std::endl;
  }

  f.close();

  if (!sfac.empty()) {
    return std::make_pair(1, "INFO:: SHELXL file " + filename + " written.");
  }
  return std::make_pair(0, "ERROR:: SHELXL file could not be written.");
}

} // namespace coot