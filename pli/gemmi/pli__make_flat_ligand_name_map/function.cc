#include "function.hh"
#include <cmath>

namespace pli {

std::map<std::string, std::string>
make_flat_ligand_name_map_gemmi(const gemmi::Residue& residue) {
    const double bond_to_H_dist = 1.1;
    const double b2Hd2 = bond_to_H_dist * bond_to_H_dist;
    std::map<std::string, std::string> map;
    
    for (const auto& atom_i : residue.atoms) {
        if (!atom_i.element.is_hydrogen()) {
            continue;
        }
        const gemmi::Position pos_i = atom_i.pos;
        
        for (const auto& atom_j : residue.atoms) {
            if (atom_i.name == atom_j.name) {
                continue;
            }
            if (atom_j.element.is_hydrogen()) {
                continue;
            }
            const gemmi::Position pos_j = atom_j.pos;
            const double dist_sq = (pos_i - pos_j).length_sq();
            if (dist_sq < b2Hd2) {
                map[atom_j.name] = atom_i.name;
                break;
            }
        }
    }
    
    return map;
}

}