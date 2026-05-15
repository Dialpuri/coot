#pragma once
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <iostream>

namespace coot {
namespace db_strands {

// Ported from: mmdb::Manager *get_mol(const std::string &filename)
// Uses gemmi::read_pdb_file instead of MMDB ReadCoorFile.
gemmi::Structure *get_mol_gemmi(const std::string &filename);

} // namespace db_strands
} // namespace coot