#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_his_ring_H__9cb67e/gemmi/function.hh"

namespace coot {
namespace reduce {

inline void
add_his_ring_H_gemmi(const std::string &H_at_name,
                     const std::string &first_neigh,
                     const std::vector<std::string> &second_neighb_vec,
                     double bl,
                     gemmi::Residue &residue) {

   if (second_neighb_vec.size() == 2) {
      coot::reduce::add_his_ring_H_gemmi(H_at_name, second_neighb_vec[0], first_neigh, second_neighb_vec[1], bl, residue);
   }
}

} // namespace reduce
} // namespace coot