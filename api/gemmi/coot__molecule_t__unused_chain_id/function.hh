#pragma once
#include <string>
#include <utility>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Return unused chain ID by checking which letters in A-Z,a-z are not used
// Takes a gemmi::Model (equivalent to mmdb::Model in original)
std::pair<bool, std::string>
unused_chain_id_gemmi(const gemmi::Model& model);

}
}