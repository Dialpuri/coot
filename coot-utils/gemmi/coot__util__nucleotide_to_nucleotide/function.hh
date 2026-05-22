#pragma once

#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <gemmi/unitcell.hpp>
#include <clipper/clipper.h>

namespace coot {
namespace util {

// Port of nucleotide_to_nucleotide from MMDB to gemmi
// Aligns nucleotide bases by matching atom positions
// Returns (success_flag, rotation+translation operator)
std::pair<bool, clipper::RTop_orth>
nucleotide_to_nucleotide_gemmi(const gemmi::Residue& residue,
                               const gemmi::Residue& std_base,
                               bool use_old_style_naming);

}
}