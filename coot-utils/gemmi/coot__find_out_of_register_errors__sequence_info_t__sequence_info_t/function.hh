#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace find_out_of_register_errors {

class sequence_info_t {
public:
   sequence_info_t(gemmi::Residue *r, int rn, const std::string &resname) :
      residue_p(r), res_no(rn), res_name(resname) {}
   gemmi::Residue *residue_p;
   int res_no;
   std::string res_name;
};

} // namespace find_out_of_register_errors
} // namespace coot