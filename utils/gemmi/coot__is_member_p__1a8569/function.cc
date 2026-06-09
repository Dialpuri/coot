#include "function.hh"

namespace coot {

bool is_member_p_gemmi(const std::vector<std::string> &v, const std::string &a) {
    for (const auto &elem : v) {
        if (elem == a) {
            return true;
        }
    }
    return false;
}

} // namespace coot
