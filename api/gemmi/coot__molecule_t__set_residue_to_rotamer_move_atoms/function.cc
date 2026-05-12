#include "function.hh"

namespace coot {
namespace molecule_t {

int set_residue_to_rotamer_move_atoms_gemmi(gemmi::Residue& res, gemmi::Residue& moving_res) {
    int i_done = 0;
    
    // Iterate through atoms in moving_res
    for (const gemmi::Atom& mov_atom : moving_res.atoms) {
        std::string atom_name_mov = mov_atom.name;
        std::string alt_loc_mov = std::string(1, mov_atom.altloc);
        
        // Look for matching atom in res
        for (gemmi::Atom& ref_atom : res.atoms) {
            std::string atom_name_ref = ref_atom.name;
            std::string alt_loc_ref = std::string(1, ref_atom.altloc);
            
            if (atom_name_mov == atom_name_ref && alt_loc_mov == alt_loc_ref) {
                // Copy coordinates
                ref_atom.pos.x = mov_atom.pos.x;
                ref_atom.pos.y = mov_atom.pos.y;
                ref_atom.pos.z = mov_atom.pos.z;
                i_done = 1;
            }
        }
    }
    
    return i_done;
}

} // namespace molecule_t
} // namespace coot