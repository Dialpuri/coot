#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <set>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue/gemmi/function.hh"

namespace coot {

class molecule_t {
public:
  static int mutate_gemmi(gemmi::Structure &st,
                          const std::string &chain_id,
                          int seqnum,
                          char ins_code,
                          const std::string &new_res_type) {
    static const std::set<std::string> valid_residues = {
      "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    if (valid_residues.find(new_res_type) == valid_residues.end()) {
      return -1;
    }

    std::string ins_str(1, ins_code);
    gemmi::Residue *res = coot::util::get_residue_gemmi(chain_id, seqnum, ins_str, st);
    if (!res) {
      return -1;
    }

    res->name = new_res_type;
    return 1;
  }
};

} // namespace coot