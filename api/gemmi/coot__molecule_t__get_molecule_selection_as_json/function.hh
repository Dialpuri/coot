#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// Parses a CID string of the form "//<ChainID>/<ResidueNumber>" and returns
// the chain ID and residue number, or false if the format is invalid.
bool parse_cid(const std::string &cid, std::string &chain_id, int &res_seqnum);

// Port of coot::molecule_t::get_molecule_selection_as_json using gemmi.
// Takes a gemmi::Structure and a CID string, returns JSON-like string
// containing atoms from the selected residue.
std::string get_molecule_selection_as_json_gemmi(const gemmi::Structure &st, const std::string &cid);

} // namespace coot
