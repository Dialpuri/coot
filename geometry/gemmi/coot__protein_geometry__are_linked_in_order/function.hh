#pragma once
#include <gemmi/model.hpp>
#include <string>
#include "mmdb2/mmdb_manager.h"

// Gemmi version of are_linked_in_order
// Returns true if 'first' and 'second' match the link's residue pairs in order
namespace coot {
namespace protein_geometry {

bool are_linked_in_order_gemmi(const gemmi::CRA& first_cra,
                               const gemmi::CRA& second_cra,
                               const mmdb::Link* link);

} // namespace protein_geometry
} // namespace coot