#pragma once

#include <string>
#include "ligand/ligand.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__molecule__read_file/gemmi/function.hh"

namespace coot {

void install_ligand_gemmi(ligand& lig, std::string ligand_pdb_filename);

} // namespace coot
