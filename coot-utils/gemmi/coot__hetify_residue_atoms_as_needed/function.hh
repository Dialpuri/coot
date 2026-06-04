#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_member_p/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__hetify_residue_atoms/gemmi/function.hh"

namespace coot {

// Inline translation of coot::util::PDB_standard_residue_types (no _gemmi port exists)
static std::vector<std::string> PDB_standard_residue_types() {
  std::vector<std::string> v;
  // Standard amino acids + MSE
  v.push_back("ALA"); v.push_back("ARG"); v.push_back("ASP");
  v.push_back("ASN"); v.push_back("CYS"); v.push_back("SER");
  v.push_back("PRO"); v.push_back("PHE"); v.push_back("GLY");
  v.push_back("GLU"); v.push_back("GLN"); v.push_back("ILE");
  v.push_back("LEU"); v.push_back("TYR"); v.push_back("TRP");
  v.push_back("HIS"); v.push_back("LYS"); v.push_back("MET");
  v.push_back("VAL"); v.push_back("THR"); v.push_back("MSE");
  // RNA single-letter and two-letter forms
  v.push_back("Td"); v.push_back("Tr"); v.push_back("T");
  v.push_back("Gd"); v.push_back("Gr"); v.push_back("G");
  v.push_back("Ad"); v.push_back("Ar"); v.push_back("A");
  // DNA
  v.push_back("DG"); v.push_back("DC"); v.push_back("DA");
  v.push_back("DU"); v.push_back("DT"); v.push_back("DI");
  // Catch-all
  v.push_back("UNK"); v.push_back("N");
  return v;
}

int hetify_residue_atoms_as_needed_gemmi(gemmi::Residue* res) {
  int r = 0;
  if (res) {
    std::string resname = res->name;
    if (!is_member_p_gemmi(PDB_standard_residue_types(), resname))
      r = hetify_residue_atoms_gemmi(*res);
  }
  return r;
}

} // namespace coot