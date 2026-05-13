#pragma once
#include <map>
#include <string>
#include <gemmi/model.hpp>

namespace pli {

std::map<std::string, std::string>
make_flat_ligand_name_map_gemmi(const gemmi::Residue& residue);

}