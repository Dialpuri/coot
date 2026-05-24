#pragma once
#include <utility>
#include <vector>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

int
cis_trans_conversion_gemmi(gemmi::Residue *res_first, gemmi::Residue *res_second,
                           gemmi::Structure &mol, gemmi::Structure &standard_residues_mol);

} // namespace util
} // namespace coot