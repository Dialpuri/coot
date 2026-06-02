#pragma once

#include "coot-utils/secondary-structure-headers.hh"
#include <gemmi/model.hpp>
#include <vector>
#include <set>

// gemmi-ported version of coot::secondary_structure_header_records::get_sheet_order
std::vector<std::vector<coot::secondary_structure_header_records::strand_relation_t>>
get_sheet_order_gemmi(
    gemmi::Structure &st,
    gemmi::Model &model,
    const std::vector<std::vector<gemmi::Residue*>> &strands_with_residues);