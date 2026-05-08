#pragma once
#include <string>
#include <vector>

namespace coot {

struct molecule_t {
  std::vector<std::string> get_chain_ids_gemmi() const;
};

}
