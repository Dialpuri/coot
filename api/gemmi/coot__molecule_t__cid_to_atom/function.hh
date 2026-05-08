#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Parse a CID string "//chain/residue_num/atom_name" and return a pointer
// to the matching atom in the structure, or nullptr if not found.
// CID format: //chain/residue_num/atom_name
// If atom_name is omitted, returns the first atom in the residue.
// If cid is empty, returns the first atom in the structure.
gemmi::Atom *cid_to_atom_gemmi(const std::string &cid, gemmi::Structure &st);

} // namespace molecule_t
} // namespace coot