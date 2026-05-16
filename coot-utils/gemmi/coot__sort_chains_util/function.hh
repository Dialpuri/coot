#pragma once

#include <string>

namespace coot {

// Comparator that sorts by chain ID string (lexicographic, ascending).
// The chain pointer field from the original mmdb::Chain* version is omitted
// because the comparison only depends on the second (string) element.
inline bool
sort_chains_util_gemmi(const std::string &a, const std::string &b) {
    return a < b;
}

}  // namespace coot