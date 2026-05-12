#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct basic_dict_restraint_t {
    std::string atom_id_1_4c() const { return atom_id_1; }
    std::string atom_id_2_4c() const { return atom_id_2; }
    std::string atom_id_1;
    std::string atom_id_2;
};

struct dictionary_residue_restraints_t {
    std::vector<basic_dict_restraint_t> bond_restraint;
};

class protein_geometry {
public:
    bool atoms_match_dictionary_bond_distance_check_gemmi(
        const gemmi::Residue& residue,
        bool check_hydrogens_too_flag,
        const dictionary_residue_restraints_t& restraints) const {
        bool status = true; // good
        
        // Skip hydrogens if check_hydrogens_too_flag is false
        auto should_include_atom = [&](const gemmi::Atom& atom) {
            if (check_hydrogens_too_flag) {
                return true;
            }
            return !atom.element.is_hydrogen();
        };
        
        // Check each bond restraint
        for (size_t ibond = 0; ibond < restraints.bond_restraint.size(); ibond++) {
            const auto& bond = restraints.bond_restraint[ibond];
            
            // Find matching atom pairs
            for (size_t iat = 0; iat < residue.atoms.size(); iat++) {
                const gemmi::Atom& at_1 = residue.atoms[iat];
                std::string atom_name_1 = at_1.name;
                
                if (bond.atom_id_1_4c() == atom_name_1 && should_include_atom(at_1)) {
                    for (size_t jat = iat + 1; jat < residue.atoms.size(); jat++) {
                        const gemmi::Atom& at_2 = residue.atoms[jat];
                        std::string atom_name_2 = at_2.name;
                        
                        if (bond.atom_id_2_4c() == atom_name_2 && should_include_atom(at_2)) {
                            // Check alt loc (in gemmi, altloc is char, not string)
                            if (at_1.altloc == at_2.altloc) {
                                // Calculate squared distance using gemmi's Position
                                double dx = at_1.pos.x - at_2.pos.x;
                                double dy = at_1.pos.y - at_2.pos.y;
                                double dz = at_1.pos.z - at_2.pos.z;
                                double d_sq = dx*dx + dy*dy + dz*dz;
                                
                                if (d_sq > 10.0) {
                                    status = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (!status) break;
            }
            
            if (!status) break;
        }
        
        return status;
    }
};

} // namespace coot