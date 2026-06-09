#pragma once

#include <vector>
#include <algorithm>

namespace coot {

template <typename T>
bool is_member_p_gemmi(const std::vector<T> &v, const T &a) {
    return std::find(v.begin(), v.end(), a) != v.end();
}

} // namespace coot
