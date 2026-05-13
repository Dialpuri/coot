#pragma once
#include <string>
#include <vector>
#include <set>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace molecule_t {

// Port of coot::molecule_t::get_residue_sidechain_average_position to gemmi
// Takes a CID string (e.g., "//A/1") and a gemmi::Model reference
// Returns a vector of 3 doubles for the average sidechain position, or empty vector if none
inline std::vector<double>
get_residue_sidechain_average_position_gemmi(const std::string &cid, gemmi::Model &model) {
   std::vector<double> v;
   gemmi::Residue *residue_p = coot::cid_to_residue_gemmi(cid, model);
   if (residue_p) {
      std::vector<clipper::Coord_orth> side_chain_positions;
      // Main chain atom names as they appear in gemmi (no trailing spaces)
      std::set<std::string> main_chain_atoms;
      main_chain_atoms.insert("CA");
      main_chain_atoms.insert("C");
      main_chain_atoms.insert("N");
      main_chain_atoms.insert("O");
      main_chain_atoms.insert("H");
      main_chain_atoms.insert("HA");
      
      for (auto &atom : residue_p->atoms) {
         std::string atom_name = atom.name;
         if (main_chain_atoms.find(atom_name) == main_chain_atoms.end()) {
            clipper::Coord_orth p = coot::co_gemmi(&atom);
            side_chain_positions.push_back(p);
         }
      }
      
      if (!side_chain_positions.empty()) {
         clipper::Coord_orth sum(0,0,0);
         for (const auto &pos : side_chain_positions)
            sum += pos;
         double is = 1.0 / static_cast<double>(side_chain_positions.size());
         v = {sum.x() * is, sum.y() * is, sum.z() * is};
      }
   }
   return v;
}

} // namespace molecule_t
} // namespace coot