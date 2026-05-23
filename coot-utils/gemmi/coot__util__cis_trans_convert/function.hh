#pragma once
#include <utility>
#include <vector>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

int
cis_trans_convert_gemmi(std::pair<gemmi::Residue*, gemmi::Residue*> mol_residues,
                        std::vector<gemmi::Residue*>& trans_residues,
                        std::vector<gemmi::Residue*>& cis_residues);

} // namespace util
} // namespace coot