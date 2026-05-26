#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

int pepflip_gemmi(gemmi::Structure &st,
                  const std::string &chain_id,
                  int resno,
                  const std::string &ins_code,
                  const std::string &altconf);

} // namespace coot