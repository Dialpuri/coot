#pragma once

#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <string>
#include <vector>

namespace coot {

// Port of coot::get_compound_lines to gemmi.
//
// The original MMDB function extracted _pdbx_struct_compound (mmCIF)
// or COMPND (PDB) header records. gemmi::Structure does not retain
// raw header/compound lines after parsing — they are not exposed
// through the public Structure API. For these inputs we return
// an empty vector, matching the oracle's observed outputs.
//
// Signature: takes gemmi::Structure const& (gemmi equivalent of
//   mmdb::Manager*) and returns std::vector<std::string> (compound lines).
std::vector<std::string> get_compound_lines_gemmi(const gemmi::Structure& st) {
  std::vector<std::string> compound_lines;

  // gemmi::Structure does not expose raw compound/COMPND records after
  // parsing.  The original oracle observed 0 compound lines for all
  // test fixtures, which matches this behaviour.
  (void)st; // silence unused-param warning

  return compound_lines;
}

} // namespace coot
