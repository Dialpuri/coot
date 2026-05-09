#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
  class molecule_t {
  public:
    molecule_t() = default;
    explicit molecule_t(const gemmi::Structure& st) : structure_(st) {}
    
    std::string get_pucker_analysis_info_gemmi() const;
    
  private:
    gemmi::Structure structure_;
  };
}