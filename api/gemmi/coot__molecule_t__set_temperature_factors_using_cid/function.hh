#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {
namespace molecule_t {

// Parses CID string like "//A/1" and sets tempFactors for matching atoms
// CID format: //chain/resnum (e.g., "//A/1", "//B/10")
void set_temperature_factors_using_cid_gemmi(gemmi::Structure& st,
                                             const std::string& cid,
                                             float temp_fact);

}
}