#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {
namespace atom_quad {

struct atom_quad {
    gemmi::Atom* atom_1;
    gemmi::Atom* atom_2;
    gemmi::Atom* atom_3;
    gemmi::Atom* atom_4;
    std::string name;
    
    atom_quad() : atom_1(nullptr), atom_2(nullptr), atom_3(nullptr), atom_4(nullptr) {}
};

// Port of coot::atom_quad::setup_chiral_quad from MMDB to gemmi
inline atom_quad setup_chiral_quad_gemmi(const gemmi::Residue* residue_with_O,
                                         const gemmi::Residue* residue_with_chiral_centre,
                                         const std::string& O_name,
                                         const std::vector<std::string>& chiral_atom_names) {
    atom_quad quad;
    
    if (chiral_atom_names.size() < 3) {
        return quad;
    }
    
    std::string atom_2_name = chiral_atom_names[0];
    std::string atom_3_name = chiral_atom_names[1];
    std::string atom_4_name = chiral_atom_names[2];
    
    if (O_name != "") { 
        for (const auto& atom : residue_with_O->atoms) {
            if (atom.name == O_name) {
                if (!quad.atom_1) {
                    quad.atom_1 = const_cast<gemmi::Atom*>(&atom);
                    break;
                }
            }
        }
        
        for (const auto& atom : residue_with_chiral_centre->atoms) {
            if (atom.name == atom_4_name) {
                if (!quad.atom_4)
                    quad.atom_4 = const_cast<gemmi::Atom*>(&atom);
            }
            if (atom.name == atom_2_name) {
                if (!quad.atom_2)
                    quad.atom_2 = const_cast<gemmi::Atom*>(&atom);
            }
            if (atom.name == atom_3_name) {
                if (!quad.atom_3)
                    quad.atom_3 = const_cast<gemmi::Atom*>(&atom);
            }
        }
    }
    return quad;
}

} // namespace atom_quad
} // namespace coot