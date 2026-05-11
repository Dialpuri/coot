#pragma once
#include <string>
#include <vector>
#include <utility>

namespace coot {

class contact_info {
public:
    std::vector<std::pair<std::string, double> > atom_radii;
    
    // Original method ported with _gemmi suffix
    double get_radius_gemmi(const std::string &element) const;
    
    // Helper to populate atom_radii for testing
    void setup_atom_radii();
};

} // namespace coot