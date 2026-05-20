#pragma once
#include <string>
#include <stdexcept>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

struct pucker_analysis_info_t {
    // Member variables that were present in the original class
    gemmi::Atom* C1_prime = nullptr;
    gemmi::Atom* N1_or_9 = nullptr;
    gemmi::Atom* O4_prime = nullptr;
    gemmi::Atom* C4_prime = nullptr;
    gemmi::Atom* C3_prime = nullptr;
    gemmi::Atom* C2_prime = nullptr;
    char altconf = ' ';
    
    // Constructor
    pucker_analysis_info_t(const gemmi::Residue& res, char altconf_char = ' ')
        : altconf(altconf_char) {
        
        for (const auto& atom : res.atoms) {
            std::string atm_name = atom.name;
            char atm_altloc = atom.altloc;
            
            // Handle alt loc - use atom if altconf is ' ' or matches atom's altloc
            if (altconf == ' ' || altconf == atm_altloc) {
                if (atm_name == " C1'" || atm_name == " C1*") {
                    C1_prime = const_cast<gemmi::Atom*>(&atom);
                } else if (atm_name == " N1 " || atm_name == " N9 ") {
                    N1_or_9 = const_cast<gemmi::Atom*>(&atom);
                } else if (atm_name == " O4'") {
                    O4_prime = const_cast<gemmi::Atom*>(&atom);
                } else if (atm_name == " C4'") {
                    C4_prime = const_cast<gemmi::Atom*>(&atom);
                } else if (atm_name == " C3'") {
                    C3_prime = const_cast<gemmi::Atom*>(&atom);
                } else if (atm_name == " C2'") {
                    C2_prime = const_cast<gemmi::Atom*>(&atom);
                }
            }
        }
        
        // Check if we found the required atoms for nucleotide
        if (!C1_prime || !N1_or_9) {
            throw std::runtime_error("Required atoms C1'/C1* and N1/N9 not found in residue");
        }
    }
    
    // Ported function
    float phosphate_distance(const gemmi::Residue& following_res) const {
        if (!C1_prime) {
            throw std::runtime_error("C1*/C1' not found in this residue");
        }
        if (!N1_or_9) {
            throw std::runtime_error("N1/N9 not found in this residue");
        }
        
        float d = 0.0f;
        bool found = false;
        
        for (const auto& atom : following_res.atoms) {
            std::string atm_name = atom.name;
            char alt_name = atom.altloc;
            
            if (atm_name == " P  ") {
                if (altconf == alt_name) {
                    gemmi::Position P_pt(atom.pos);
                    gemmi::Position N_pt(N1_or_9->pos);
                    gemmi::Position C_pt(C1_prime->pos);
                    
                    gemmi::Vec3 CN = N_pt - C_pt;
                    gemmi::Vec3 CP = P_pt - C_pt;
                    
                    double CN_d = sqrt(CN.length_sq());
                    double CP_d = sqrt(CP.length_sq());
                    
                    if (CN_d > 0.0) {
                        if (CP_d > 0.0) {
                            found = true;
                            double cos_alpha = CN.dot(CP) / (CN_d * CP_d);
                            double alpha = acos(cos_alpha);
                            double sin_pi_minus_alpha = sin(M_PI - alpha);
                            d = sin_pi_minus_alpha * static_cast<float>(CP_d);
                        }
                    }
                }
            }
        }
        
        if (!found) {
            throw std::runtime_error("P not found in this residue");
        }
        
        return d;
    }
};

} // namespace coot