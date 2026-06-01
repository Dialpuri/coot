#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {

// Forward declarations of types from coot that we reuse
struct atom_name_pair {
  std::string atom1;
  std::string atom2;
  atom_name_pair() = default;
  atom_name_pair(const std::string& a1, const std::string& a2) : atom1(a1), atom2(a2) {}
};

struct atom_index_pair {
  int index1;
  int index2;
  atom_index_pair() : index1(-1), index2(-1) {}
  atom_index_pair(int i1, int i2) : index1(i1), index2(i2) {}
};

namespace monomer_utils {

std::vector<atom_index_pair>
get_atom_index_pairs_gemmi(
  const std::vector<atom_name_pair> &atom_name_pairs_in,
  const gemmi::Residue &residue);

}  // namespace monomer_utils
}  // namespace coot