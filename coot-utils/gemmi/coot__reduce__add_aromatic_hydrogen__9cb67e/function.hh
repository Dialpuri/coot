#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"

namespace coot { namespace reduce {

inline void
add_aromatic_hydrogen_gemmi(const std::string &H_at_name,
                            const std::string &neighb_at_name_1,
                            const std::string &neighb_at_name_2,
                            const std::string &neighb_at_name_3,
                            double bl, gemmi::Residue& residue) {
   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
   for (const auto& alt_conf : alt_confs) {
      gemmi::Atom* at_n_1 = nullptr;
      gemmi::Atom* at_n_2 = nullptr;
      gemmi::Atom* at_n_3 = nullptr;
      
      for (auto& atom : residue.atoms) {
         char alt = atom.altloc;
         std::string atom_alt_conf = (alt == ' ' || alt == '\0') ? "" : std::string(1, alt);
         if (atom_alt_conf == alt_conf) {
            if (atom.name == neighb_at_name_1) at_n_1 = &atom;
            if (atom.name == neighb_at_name_2) at_n_2 = &atom;
            if (atom.name == neighb_at_name_3) at_n_3 = &atom;
         }
      }
      
      if (at_n_1 && at_n_2 && at_n_3) {
         float bf = at_n_2->b_iso;
         clipper::Coord_orth H_pos = position_by_bisection_gemmi(at_n_1, at_n_2, at_n_3, bl);
         add_hydrogen_atom_gemmi(H_at_name, H_pos, bf, alt_conf, residue);
      }
   }
}

}} // namespace coot::reduce