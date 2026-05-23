#pragma once
#include <string>
#include <iostream>

namespace coot {

struct molecule_t {
    // dummy — actual definition in coot source; we only need the method signature
    void print_secondary_structure_info() const;
    // the gemmi variant
    void print_secondary_structure_info_gemmi() const;
};

namespace util {

// print secondary structure info to std::cout (same behavior as original)
void print_secondary_structure_info(void *model_ptr);

} // namespace util
} // namespace coot