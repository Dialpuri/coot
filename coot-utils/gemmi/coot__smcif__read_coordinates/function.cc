#include "function.hh"
#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include <gemmi/elem.hpp>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>

namespace coot::smcif {

static std::pair<std::string, int> symbol_to_element(const std::string& symbol);

std::vector<gemmi::Atom> read_coordinates_gemmi(
    const std::string& cif_path,
    const clipper::Cell& cell,
    const clipper::Spacegroup& spg) {

  std::vector<gemmi::Atom> atom_vec;

  gemmi::cif::Document doc = gemmi::read_cif_gz(cif_path);
  if (doc.blocks.empty()) {
    return atom_vec;
  }

  gemmi::cif::Block& bloc = doc.blocks[0];

  // Find _atom_site loop using the label column
  gemmi::cif::Column label_col = bloc.find_loop("_atom_site_label");
  if (label_col.length() <= 0) {
    return atom_vec;
  }

  // Get other columns from the same loop
  gemmi::cif::Column fract_x_col = bloc.find_loop("_atom_site_fract_x");
  gemmi::cif::Column fract_y_col = bloc.find_loop("_atom_site_fract_y");
  gemmi::cif::Column fract_z_col = bloc.find_loop("_atom_site_fract_z");

  if (fract_x_col.length() <= 0 || fract_y_col.length() <= 0 || fract_z_col.length() <= 0) {
    return atom_vec;
  }

  // Optional columns
  gemmi::cif::Column type_symbol_col = bloc.find_loop("_atom_site_type_symbol");
  gemmi::cif::Column occupancy_col = bloc.find_loop("_atom_site_occupancy");
  gemmi::cif::Column u_iso_col = bloc.find_loop("_atom_site_U_iso_or_equiv");
  gemmi::cif::Column symm_mult_col = bloc.find_loop("_atom_site_symmetry_multiplicity");
  gemmi::cif::Column disorder_col = bloc.find_loop("_atom_site_pdbx_PDB_disorder_id");

  int n_rows = label_col.length();
  for (int row = 0; row < n_rows; ++row) {
    std::string label = label_col.at(row);
    double xf, yf, zf;
    try {
      xf = std::stod(fract_x_col.at(row));
      yf = std::stod(fract_y_col.at(row));
      zf = std::stod(fract_z_col.at(row));
    } catch (...) {
      continue;
    }

    // Get element symbol
    std::string symbol = "UNK";
    if (type_symbol_col.length() > 0) {
      std::string sym = type_symbol_col.at(row);
      if (!sym.empty() && sym != ".") {
        symbol = sym;
      }
    }
    if (symbol == "UNK") {
      std::string s = label;
      s.erase(std::remove_if(s.begin(), s.end(), ::isdigit), s.end());
      symbol = s.empty() ? "UNK" : s;
    }

    // Get disorder/alt loc
    std::string alt_loc;
    if (disorder_col.length() > 0) {
      alt_loc = disorder_col.at(row);
    }

    // Occupancy and B-factor
    double occ = 1.0;
    double tf = 10.0;
    int symmetry_multiplicity = 1;
    bool has_symm_mult = false;

    if (occupancy_col.length() > 0) {
      try { occ = std::stod(occupancy_col.at(row)); } catch (...) {}
    }

    if (u_iso_col.length() > 0) {
      try {
        double u_iso = std::stod(u_iso_col.at(row));
        tf = u_iso * 8.0 * M_PI * M_PI;
      } catch (...) {}
    }

    if (symm_mult_col.length() > 0) {
      try {
        symmetry_multiplicity = std::stoi(symm_mult_col.at(row));
        has_symm_mult = true;
      } catch (...) {}
    }

    // Convert fractional to orthogonal
    clipper::Coord_frac cf(xf, yf, zf);
    clipper::Coord_orth co = cf.coord_orth(cell);

    double occ_symm = occ;
    if (has_symm_mult) {
      occ_symm /= static_cast<float>(symmetry_multiplicity);
    }

    // Create atom
    gemmi::Atom at;
    at.name = label;
    at.pos = gemmi::Position(co.x(), co.y(), co.z());
    at.occ = occ_symm;
    at.b_iso = tf;

    std::pair<std::string, int> ele = symbol_to_element(symbol);
    at.element = gemmi::Element(ele.first);

    if (!alt_loc.empty()) {
      at.altloc = alt_loc[0];
    }

    atom_vec.push_back(at);
  }

  return atom_vec;
}

static std::pair<std::string, int> symbol_to_element(const std::string& symbol) {
  if (symbol == "C")  return {"C", 0};
  if (symbol == "N")  return {"N", 0};
  if (symbol == "O")  return {"O", 0};
  if (symbol == "H")  return {"H", 0};
  if (symbol == "S")  return {"S", 0};
  if (symbol == "P")  return {"P", 0};
  if (symbol == "F")  return {"F", 0};
  if (symbol == "CL") return {"Cl", 0};
  if (symbol == "CA") return {"Ca", 0};
  if (symbol == "NA") return {"Na", 0};
  if (symbol == "MG") return {"Mg", 0};
  if (symbol == "MG+2") return {"Mg", 2};
  if (symbol == "ZN") return {"Zn", 0};
  if (symbol == "FE") return {"Fe", 0};
  if (symbol == "CO") return {"Co", 0};
  if (symbol == "MN") return {"Mn", 0};
  if (symbol == "CU") return {"Cu", 0};
  if (symbol == "NI") return {"Ni", 0};
  if (symbol == "LI") return {"Li", 0};
  if (symbol == "K")  return {"K", 0};
  if (symbol == "RU") return {"Ru", 0};
  if (symbol == "CS") return {"Cs", 0};
  if (symbol == "SR") return {"Sr", 0};
  if (symbol == "BR") return {"Br", 0};
  if (symbol == "I")  return {"I", 0};
  if (symbol == "SE") return {"Se", 0};
  std::string elem;
  if (!symbol.empty()) {
    elem += static_cast<char>(std::toupper(symbol[0]));
    if (symbol.size() > 1 && std::isalpha(static_cast<unsigned char>(symbol[1]))) {
      elem += static_cast<char>(std::tolower(symbol[1]));
    }
    return {elem, 0};
  }
  return {"?", 0};
}

} // namespace coot::smcif