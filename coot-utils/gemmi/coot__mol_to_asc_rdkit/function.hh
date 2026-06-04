#pragma once

#include <gemmi/model.hpp>
#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <unordered_map>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__make_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__create_mmdbmanager_from_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/make_asc/gemmi/function.hh"

namespace coot {

atom_selection_container_t mol_to_asc_rdkit_gemmi(const std::string& file_name);

} // namespace coot