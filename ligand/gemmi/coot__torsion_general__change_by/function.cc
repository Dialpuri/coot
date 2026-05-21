#include "function.hh"
#include <cmath>
#include <gemmi/math.hpp>

namespace coot {
namespace torsion_general {

int change_by_gemmi(double diff_degrees,
                    std::vector<gemmi::Position>& coords,
                    int atom1_idx, int atom2_idx,
                    const std::vector<int>& children_of_atom2) {
    using namespace gemmi;
    
    // Convert degrees to radians
    double diff_rad = diff_degrees * M_PI / 180.0;
    
    // Get the bond axis vector (from atom1 to atom2)
    Vec3 bond_axis = coords[atom2_idx] - coords[atom1_idx];
    double bond_len = bond_axis.length();
    if (bond_len < 1e-6)
        return 1;  // degenerate case
    
    // Normalize the bond axis
    Vec3 axis = bond_axis / bond_len;
    
    // Transform all atoms that are children of atom2
    // (in the original code, these are atoms that will rotate when we rotate about the bond)
    // Note: atom2 is the second atom in the bond (clicked_atom_indices[2] in original)
    // The children are the atoms that are descendants of atom2 in the tree
    for (size_t i = 0; i < coords.size(); ++i) {
        bool should_rotate = false;
        
        // Check if this atom is in the subtree rooted at atom2
        for (int child_idx : children_of_atom2) {
            if (static_cast<int>(i) == child_idx) {
                should_rotate = true;
                break;
            }
        }
        
        if (should_rotate) {
            // Transform the position: rotate around atom2's position using rotate_about_axis
            Vec3 relative_pos = static_cast<Vec3>(coords[i] - coords[atom2_idx]);
            Vec3 rotated_pos = rotate_about_axis(relative_pos, axis, diff_rad);
            coords[i] = Position(static_cast<Vec3>(coords[atom2_idx]) + rotated_pos);
        }
    }
    
    return 0;  // success
}

} // namespace torsion_general
} // namespace coot