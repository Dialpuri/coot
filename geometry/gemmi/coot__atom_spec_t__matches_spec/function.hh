#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {

class atom_spec_t {
public:
  std::string chain_id;
  int res_no;
  std::string ins_code;
  std::string atom_name;
  std::string alt_conf;
  int int_user_data;
  float float_user_data;
  std::string string_user_data;
  int model_number;

  atom_spec_t() : chain_id("unset") {
    res_no = -100000;
    model_number = -1;
    int_user_data = -1;
    float_user_data = -1;
  }

  bool matches_spec_gemmi(gemmi::CRA cra) const {
    if (cra.atom == nullptr)
      return false;

    if (atom_name == cra.atom->name) {
      if (alt_conf == std::string(1, cra.atom->altloc)) {
        if (cra.residue) {
          if (res_no == cra.residue->seqid.num.value) {
            char icode = cra.residue->seqid.icode;
            std::string res_icode = std::string(1, icode);
            if (ins_code == res_icode) {
              if (cra.chain) {
                if (chain_id == cra.chain->name) {
                  return true;
                } else {
                  return false;
                }
              } else {
                return true;
              }
            } else {
              return false;
            }
          } else {
            return false;
          }
        } else {
          return true;
        }
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
};

}