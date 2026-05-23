#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

struct atom_spec_t {
  std::string chain_id;
  int res_no;
  std::string ins_code;
  std::string atom_name;
  std::string alt_conf;
  int model_number;
  
  atom_spec_t() : chain_id(""), res_no(0), ins_code(""), atom_name(""), alt_conf(""), model_number(0) {}
  atom_spec_t(const std::string& chain_in, int resno_in, const std::string& insertion_code_in,
              const std::string& atom_name_in, const std::string& alt_conf_in)
    : chain_id(chain_in), res_no(resno_in), ins_code(insertion_code_in),
      atom_name(atom_name_in), alt_conf(alt_conf_in), model_number(0) {}
};

namespace torsion_general {

// atom_index_gemmi - find index of atom matching spec in residue
// Returns -1 (FAIL) if not found, otherwise 0-based index into residue.atoms
int atom_index_gemmi(const atom_spec_t& spec, gemmi::Residue& residue);

} // namespace torsion_general
} // namespace coot