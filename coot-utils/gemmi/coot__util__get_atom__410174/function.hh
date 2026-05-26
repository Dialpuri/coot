#pragma once
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue__ffd79d/gemmi/function.hh"

namespace coot { namespace util {

// Get atom from gemmi Structure based on atom spec parameters
inline gemmi::Atom *
get_atom_gemmi(const std::string &chain_id, int res_no, const std::string &ins_code,
               const std::string &atom_name, const std::string &alt_conf,
               gemmi::Structure &st) {
   gemmi::Residue *res = get_residue_gemmi(chain_id, res_no, ins_code, st);
   if (!res) return nullptr;

   for (gemmi::Atom &atom : res->atoms) {
      if (atom.name == atom_name) {
         // Handle altloc comparison - gemmi uses altloc char, empty string in spec means any
         if (alt_conf.empty() || atom.altloc == alt_conf[0]) {
            return &atom;
         }
      }
   }
   return nullptr;
}

}}