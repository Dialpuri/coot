#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"

namespace coot { namespace reduce {

namespace {

static std::string trim_ws(const std::string& s) {
    auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

} // anonymous namespace

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

      std::string trimmed_1 = trim_ws(neighb_at_name_1);
      std::string trimmed_2 = trim_ws(neighb_at_name_2);
      std::string trimmed_3 = trim_ws(neighb_at_name_3);

      for (auto& atom : residue.atoms) {
         char alt = atom.altloc;
         std::string atom_alt_conf = (alt == ' ' || alt == '\0') ? "" : std::string(1, alt);
         if (atom_alt_conf == alt_conf) {
            std::string atom_name_trimmed = trim_ws(atom.name);
            if (atom_name_trimmed == trimmed_1) at_n_1 = &atom;
            if (atom_name_trimmed == trimmed_2) at_n_2 = &atom;
            if (atom_name_trimmed == trimmed_3) at_n_3 = &atom;
         }
      }

      if (at_n_1 && at_n_2 && at_n_3) {
         float bf = at_n_2->b_iso;
         clipper::Coord_orth H_pos = position_by_bisection_gemmi(at_n_1, at_n_2, at_n_3, bl);
         add_hydrogen_atom_gemmi(H_at_name, H_pos, bf, alt_conf, residue);
      }
   }
}

inline void
add_aromatic_hydrogen_gemmi(const std::string &H_at_name,
                            const std::string &first_neigh,
                            const std::vector<std::string> &second_neighb_vec,
                            double bl,
                            gemmi::Residue *residue_p) {
   if (second_neighb_vec.size() == 2) {
      add_aromatic_hydrogen_gemmi(H_at_name, second_neighb_vec[0], first_neigh, second_neighb_vec[1],
                                  bl, *residue_p);
   }
}

}} // namespace coot::reduce
