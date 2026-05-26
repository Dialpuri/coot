#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class secondary_structure_header_records {
public:
    secondary_structure_header_records() = default;
    secondary_structure_header_records(const secondary_structure_header_records& other) = default;
    secondary_structure_header_records& operator=(const secondary_structure_header_records& other) = default;
    
    // Constructor that takes a gemmi Structure instead of MMDB manager
    secondary_structure_header_records(const gemmi::Structure& /*st*/, bool /*flag*/) {}
};

} // namespace coot