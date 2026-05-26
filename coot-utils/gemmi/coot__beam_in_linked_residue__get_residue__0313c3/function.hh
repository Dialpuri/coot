#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

namespace beam_in_linked_residue {

// Port of mmdb::Residue* get_residue(const std::string &comp_id, mmdb::Manager *mol) const
// Returns gemmi::Residue* instead of mmdb::Residue*
gemmi::Residue* get_residue_gemmi(const std::string &comp_id, const gemmi::Structure &st);

} // namespace beam_in_linked_residue
} // namespace coot