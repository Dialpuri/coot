#include "function.hh"
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <iostream>
#include <gemmi/elem.hpp>

namespace coot {

float ShelxIns::string_to_float(const std::string& s) {
    return std::stof(s);
}

std::string ShelxIns::make_atom_element(const std::string& name, int sfac_index) const {
    // sfac_index is 1-based in shelx (from SFAC table)
    if (sfac_index <= 0 || sfac_index > static_cast<int>(sfac.size())) {
        return "ERROR-in-SFAC";
    }
    return sfac[sfac_index - 1];
}

std::string ShelxIns::make_atom_name(const std::string& name, const std::string& element) const {
    // Format the atom name to match MMDB's GetAtomName format (4 chars)
    // MMDB pads with space at beginning if name is 3 chars or less
    // Example: "C1" becomes " C1 " (1 leading space, 1 trailing space)
    if (name.length() == 4) {
        return name;
    } else if (name.length() < 4) {
        std::string result = " " + name;
        while (result.length() < 4) {
            result += " ";
        }
        return result;
    } else {
        return name.substr(0, 4);
    }
}

gemmi::Atom* ShelxIns::previous_non_riding_atom(const std::vector<gemmi::Atom*>& atom_vector,
                                               int udd_non_riding_atom_flag_handle_in) const {
    for (auto it = atom_vector.rbegin(); it != atom_vector.rend(); ++it) {
        gemmi::Atom* atom = *it;
        // Check if this atom is marked as non-riding
        // In the original code, this uses PutUDData to store flag value
        // For gemmi port, we'll check a simple heuristic: if tempFactor > 0, it's non-riding
        if (atom->b_iso > 0) {
            return atom;
        }
    }
    return nullptr;
}

static std::string format_element_for_mmdb(const std::string& elem) {
    // MMDB's GetElementName returns 2-character padded string
    // Example: "C" becomes " C" (1 leading space, 1 char for element)
    if (elem.length() == 1) {
        return " " + elem;
    } else if (elem.length() == 2) {
        return elem;
    } else {
        return elem.substr(0, 2);
    }
}

gemmi::Atom* ShelxIns::make_atom_gemmi(const shelx_card_info_t& card, const std::string& altconf,
                                       int udd_afix_handle_in,
                                       int udd_non_riding_atom_flag_handle_in,
                                       int udd_riding_atom_negative_u_value_handle_in,
                                       bool have_udd_atoms, int current_afix,
                                       const clipper::Cell& cell_in,
                                       const std::vector<gemmi::Atom*>& atom_vector) const {
    
    gemmi::Atom* at = new gemmi::Atom;  // returned
    
    // local
    int sfac_index;
    double u_to_b = 8.0 * M_PI * M_PI;  // perhaps this should be a function
    
    // std::cout << "DEBUG:: new atom for " <<  card.card << std::endl;
    
    sfac_index = std::atoi(card.words[1].c_str());
    std::string element = make_atom_element(card.words[0], sfac_index);
    if (element == "ERROR-in-SFAC") {
        std::cout << "WARNING:: problem making element - ignoring atom" << std::endl;
        delete at;
        at = nullptr;
    } else { 
        at->name = make_atom_name(card.words[0], element);
        at->pos.x = std::atof(card.words[2].c_str());
        at->pos.y = std::atof(card.words[3].c_str());
        at->pos.z = std::atof(card.words[4].c_str());
        float occupancy = 1.0;
        float b_synth = 10.0;

        try {
            if (card.words.size() > 5)
                occupancy = string_to_float(card.words[5]);
            
            at->pos.x = string_to_float(card.words[2]);
            at->pos.y = string_to_float(card.words[3]);
            at->pos.z = string_to_float(card.words[4]);
            at->occ = occupancy;
            at->b_iso = b_synth;
            at->element = gemmi::Element(element.c_str());
            if (!altconf.empty()) {
                at->altloc = altconf[0];
            }
        }
        catch (const std::runtime_error &rte) {
            // do nothing
        }
        if (card.words.size() >= 6) {
            
            // What makes an atom non-riding?
            // 
            // I'm guessing that there is an anisotropic U or sing U > 0
            
            if (card.words.size() > 6) {
                if (card.words.size() < 8) {
                    // isotropic temperature factor
                    double u_factor_from_card = std::atof(card.words[6].c_str());
                    if (u_factor_from_card > 0.0 ) { 
                        at->b_iso = u_to_b * u_factor_from_card;
                        at->flag = at->flag | 1;  // basic flag
                        // For gemmi, we don't have UDD handles, so skip that
                    } else {
                        // negative U:
                        // 
                        // If U is between -0.5 and -5.0 then apply riding U rule
                        // 
                        if ((u_factor_from_card <= -0.5) && (u_factor_from_card >= -5.0)) {
                            // Find previous non-riding atom and use that to determine b for this atom

                            gemmi::Atom *prev = previous_non_riding_atom(atom_vector, udd_non_riding_atom_flag_handle_in);
                            if (prev) {
                                // Store the negative U value for reference
                                // For gemmi, we could use a custom property
                                at->b_iso = prev->b_iso * -u_factor_from_card;
                            } else {
                                // Don't know what to do.  Does this ever happen?
                                at->b_iso = u_factor_from_card;
                            } 
                        } else {
                            // Don't know what to do.  Does this ever happen?
                            at->b_iso = u_factor_from_card;
                        }
                    }
                } else {
                    if (card.words.size() > 11) {
                        // anisotropic temperature factor
                        double u11 = std::atof(card.words[6].c_str());
                        double u22 = std::atof(card.words[7].c_str());
                        double u33 = std::atof(card.words[8].c_str());
                        double u23 = std::atof(card.words[9].c_str());
                        double u13 = std::atof(card.words[10].c_str());
                        double u12 = std::atof(card.words[11].c_str());

                        double a = cell_in.a();
                        double b = cell_in.b();
                        double c = cell_in.c();

                        // Now othogonalize the U values:
                        // clipper::U_aniso_frac ocaf(u11, u22, u33, u12, u13, u23);
                        clipper::U_aniso_frac caf(u11/(a*a), u22/(b*b), u33/(c*c),
                                                  u12/(a*b), u13/(a*c), u23/(b*c));
                        clipper::U_aniso_orth cao = caf.u_aniso_orth(cell_in);

                        // Store in gemmi aniso format: u11, u22, u33, u23, u13, u12
                        at->aniso.u11 = cao(0,0);
                        at->aniso.u22 = cao(1,1);
                        at->aniso.u33 = cao(2,2);
                        at->aniso.u23 = cao(1,2);
                        at->aniso.u13 = cao(0,2);
                        at->aniso.u12 = cao(0,1);

                        at->flag = at->flag | 1;  // basic flag
                        double u_synth = (at->aniso.u11 + at->aniso.u22 + at->aniso.u33)/3.0;
                        at->b_iso = u_to_b * u_synth;
                        // For gemmi, we don't have UDD handles
                    }
                }
            } else {
                // An atom with minimal description. Let's make up a
                // temperature factor:
                at->b_iso = 1.0;
                at->flag = at->flag | 1;  // basic flag
            } 
        
            // For gemmi, we don't have UDD handles
        } else {
            std::cout << "(make_atom) bad atom: " << card.card << std::endl;
            delete at;
            at = nullptr;
        }
    }
    return at;
}

} // namespace coot