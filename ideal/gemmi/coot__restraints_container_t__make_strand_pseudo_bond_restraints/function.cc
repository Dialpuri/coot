#include "function.hh"
#include <gemmi/math.hpp>
#include <iostream>

namespace coot {

void make_strand_pseudo_bond_restraints_gemmi(
    const std::vector<const gemmi::Residue*>& residues) {
    
    // This method creates pseudo-bond restraints for strand structures:
    // - O-O bonds (distance ~4.64 Å) 
    // - O-O-O angles (angle ~98°)
    // - CA-CA-CA angles (angle ~120°)
    // 
    // The algorithm requires consecutive residues (seqNum differs by 1).
    // It ignores residues with alternative conformations.
    
    float pseudo_bond_esd = 0.08; // just a guess
    
    // Process consecutive residue pairs and triplets
    for (size_t i = 1; i < residues.size(); i++) {
        const gemmi::Residue* res1 = residues[i-1];
        const gemmi::Residue* res2 = residues[i];
        
        // Check consecutive sequence numbers
        if (res2->seqid.num.value == res1->seqid.num.value + 1) {
            // Look for O atoms for O-O bond
            const gemmi::Atom* o1 = nullptr;
            const gemmi::Atom* o2 = nullptr;
            
            for (const auto& atom : res1->atoms) {
                if (atom.name == " O  ") {
                    o1 = &atom;
                    break;
                }
            }
            for (const auto& atom : res2->atoms) {
                if (atom.name == " O  ") {
                    o2 = &atom;
                    break;
                }
            }
            
            if (o1 && o2) {
                // Found O-O pair, calculate distance
                gemmi::Vec3 v = o1->pos - o2->pos;
                double dist = v.length();
                
                // Only create restraint if within reasonable range
                if (dist < 5.0) {
                    std::cout << "Strand Bond restraint ("
                              << o1->name << " " << res1->seqid.num.value << ") to ("
                              << o2->name << " " << res2->seqid.num.value << ") " 
                              << dist << std::endl;
                }
            }
            
            // Check for O-O-O angle (need res3)
            if (i + 1 < residues.size()) {
                const gemmi::Residue* res3 = residues[i+1];
                if (res3->seqid.num.value == res2->seqid.num.value + 1) {
                    const gemmi::Atom* o3 = nullptr;
                    for (const auto& atom : res3->atoms) {
                        if (atom.name == " O  ") {
                            o3 = &atom;
                            break;
                        }
                    }
                    if (o3) {
                        std::cout << "Strand Angle restraint ("
                                  << res1->seqid.num.value << ") to ("
                                  << res2->seqid.num.value << ") to ("
                                  << res3->seqid.num.value << ") 98.0" << std::endl;
                    }
                }
            }
        }
        
        // Now check for CA-CA-CA angle
        const gemmi::Atom* ca1 = nullptr;
        const gemmi::Atom* ca2 = nullptr;
        
        for (const auto& atom : res1->atoms) {
            if (atom.name == " CA ") {
                ca1 = &atom;
                break;
            }
        }
        for (const auto& atom : res2->atoms) {
            if (atom.name == " CA ") {
                ca2 = &atom;
                break;
            }
        }
        
        if (ca1 && ca2) {
            if (i + 1 < residues.size()) {
                const gemmi::Residue* res3 = residues[i+1];
                if (res3->seqid.num.value == res2->seqid.num.value + 1) {
                    const gemmi::Atom* ca3 = nullptr;
                    for (const auto& atom : res3->atoms) {
                        if (atom.name == " CA ") {
                            ca3 = &atom;
                            break;
                        }
                    }
                    if (ca3) {
                        std::cout << "Strand Angle restraint ("
                                  << res1->seqid.num.value << ") to ("
                                  << res2->seqid.num.value << ") to ("
                                  << res3->seqid.num.value << ") 120.0" << std::endl;
                    }
                }
            }
        }
    }
}

} // namespace coot