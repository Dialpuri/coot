#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

std::vector<std::vector<const gemmi::Chain*>>
ncs_related_chains_gemmi(const gemmi::Model& model);

std::vector<std::vector<std::string>>
get_ncs_related_chains_gemmi(const gemmi::Structure& st);

}
}