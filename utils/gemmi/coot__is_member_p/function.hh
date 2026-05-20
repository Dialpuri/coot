#pragma once
#include <vector>
#include <string>

namespace coot {

// Port of coot::is_member_p — check if string 'a' is in vector 'v'
bool is_member_p_gemmi(const std::vector<std::string> &v, const std::string &a) {
   bool ir = false;
   for (const auto& elem : v) {
      if (elem == a) {
         ir = true;
         break;
      }
   }
   return ir;
}

} // namespace coot