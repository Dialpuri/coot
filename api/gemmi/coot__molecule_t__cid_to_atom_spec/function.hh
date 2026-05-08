#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

struct atom_spec_t {
   std::string chain_id;
   int res_no;
   std::string ins_code;
   std::string atom_name;
   std::string alt_conf;
   int model_number;
   
   atom_spec_t() : chain_id("unset"), res_no(0), model_number(1) {}
};

std::pair<bool, atom_spec_t>
cid_to_atom_spec_gemmi(const gemmi::Structure& st, const std::string& cid);

}