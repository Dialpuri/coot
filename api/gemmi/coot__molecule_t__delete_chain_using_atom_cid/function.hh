#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

int
delete_chain_using_atom_cid_gemmi(gemmi::Structure& st, const std::string& cid);

}