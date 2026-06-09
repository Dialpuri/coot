#pragma once

#include <string>
#include <vector>

namespace coot {

/// Return true if string `a` is found in vector `v`.
bool is_member_p_gemmi(const std::vector<std::string> &v, const std::string &a);

} // namespace coot
