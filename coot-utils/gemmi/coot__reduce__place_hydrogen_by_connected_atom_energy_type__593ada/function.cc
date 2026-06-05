#include "function.hh"
#include <algorithm>

namespace coot::reduce {

std::vector<std::string>
place_hydrogen_by_connected_atom_energy_type_gemmi(
    const std::string& energy_type,
    unsigned int iat,
    unsigned int iat_neighb,
    const coot::dictionary_residue_restraints_t& rest,
    gemmi::Residue* residue_p,
    double bl_aliph,
    double bl_arom,
    double bl_amino,
    double bl_oh,
    double bl_sh) {

    std::vector<std::string> placed_Hs;
    if (!residue_p) return placed_Hs;

    const std::string& H_at_name  = rest.atom_info[iat].atom_id_4c;
    const std::string& at_name_1  = rest.atom_info[iat_neighb].atom_id_4c;

    double bl = bl_aliph;

    if (energy_type == "OH" || energy_type == "OS") {
        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::string at_name_2 = "";
        std::string at_name_3 = "";
        for (unsigned int i2 : iat_neighb_2) {
            if (i2 != iat) {
                at_name_2 = rest.atom_info[i2].atom_id_4c;
                break;
            }
        }
        if (!at_name_2.empty()) {
            std::vector<unsigned int> iat_neighb_3 = rest.neighbours(iat_neighb_2[0], false);
            for (unsigned int i3 : iat_neighb_3) {
                if (i3 != iat_neighb) {
                    at_name_3 = rest.atom_info[i3].atom_id_4c;
                    break;
                }
            }
        }
        bl = bl_oh;
        if (energy_type == "OS") bl = bl_sh;

        std::string other_H = get_other_H_name(iat, rest);
        bool found_other_H  = !other_H.empty();
        double tor_inital   = 0.0;
        double angle        = 109.47;

        if (found_other_H) {
            coot::reduce::add_OH_H_gemmi(
                H_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(H_at_name);

            std::string other_at_name = other_H;
            coot::reduce::add_OH_H_gemmi(
                other_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(other_at_name);
        } else {
            coot::reduce::add_OH_H_gemmi(
                H_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(H_at_name);
        }

    } else if (energy_type == "OS2") {
        std::string other_H = get_other_H_name(iat, rest);
        bool found_other_H = !other_H.empty();

        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::string at_name_2 = "";
        if (!iat_neighb_2.empty()) {
            at_name_2 = rest.atom_info[iat_neighb_2[0]].atom_id_4c;
        }
        std::vector<unsigned int> iat_neighb_3_vec = rest.neighbours(iat_neighb_2[0], false);
        std::string at_name_3 = "";
        for (unsigned int i3 : iat_neighb_3_vec) {
            if (i3 != iat_neighb) {
                at_name_3 = rest.atom_info[i3].atom_id_4c;
                break;
            }
        }

        double tor_inital = 0.0;
        double angle = 104.5;
        bl = bl_oh;

        if (found_other_H) {
            std::string other_at_name = other_H;
            coot::reduce::add_OH_H_gemmi(
                H_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(H_at_name);

            coot::reduce::add_OH_H_gemmi(
                other_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(other_at_name);
        } else {
            coot::reduce::add_OH_H_gemmi(
                H_at_name, at_name_1, at_name_2, at_name_3,
                bl, angle, tor_inital, *residue_p);
            placed_Hs.push_back(H_at_name);
        }

    } else if (energy_type == "NT3" || energy_type == "NS3" ||
               energy_type == "NC3" || energy_type == "NO2" || energy_type == "NS2") {
        // Amino nitrogen types: add ALL hydrogens connected to this nitrogen
        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::vector<std::string> second_neighb_vec;
        for (unsigned int i2 : iat_neighb_2) {
            if (i2 != iat) {
                second_neighb_vec.push_back(rest.atom_info[i2].atom_id_4c);
            }
        }

        bl = bl_amino;
        std::string at_name_2 = second_neighb_vec.size() > 0 ? second_neighb_vec[0] : "";
        std::string at_name_3 = second_neighb_vec.size() > 1 ? second_neighb_vec[1] : "";

        // Find ALL hydrogen atoms connected to iat_neighb (the nitrogen)
        // by iterating over ALL atoms in the dictionary
        std::vector<std::string> all_H_names;
        for (unsigned int k = 0; k < rest.atom_info.size(); k++) {
            if (!rest.atom_info[k].is_hydrogen()) continue;
            // Check if atom k is connected to iat_neighb
            std::vector<unsigned int> k_neighbs = rest.neighbours(k, false);
            for (unsigned int nk : k_neighbs) {
                if (nk == iat_neighb) {
                    all_H_names.push_back(rest.atom_info[k].atom_id_4c);
                    break;
                }
            }
        }

        // Put the requested H first, then the others
        std::string requested_H = rest.atom_info[iat].atom_id_4c;
        std::vector<std::string> ordered_Hs;
        ordered_Hs.push_back(requested_H);
        for (const auto& h_name : all_H_names) {
            if (h_name != requested_H) {
                ordered_Hs.push_back(h_name);
            }
        }

        // Add each H in order
        for (const auto& h_name : ordered_Hs) {
            coot::reduce::add_amino_single_H_gemmi(
                h_name, at_name_1, at_name_2, at_name_3, bl, *residue_p);
            placed_Hs.push_back(h_name);
        }

    } else if (energy_type == "C3" || energy_type == "S3") {
        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::string at_name_2 = "";
        std::string at_name_3 = "";
        int idx = 0;
        for (unsigned int i2 : iat_neighb_2) {
            if (i2 != iat) {
                if (idx == 0) at_name_2 = rest.atom_info[i2].atom_id_4c;
                else if (idx == 1) at_name_3 = rest.atom_info[i2].atom_id_4c;
                idx++;
            }
        }
        coot::reduce::add_tetrahedral_hydrogen_gemmi(
            H_at_name, at_name_1, at_name_2, at_name_3, "", bl_aliph, *residue_p);
        placed_Hs.push_back(H_at_name);

    } else if (energy_type == "C2" || energy_type == "S2") {
        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::vector<std::string> second_neighb_vec;
        for (unsigned int i2 : iat_neighb_2) {
            if (i2 != iat) {
                second_neighb_vec.push_back(rest.atom_info[i2].atom_id_4c);
            }
        }

        std::string H_at_name_2 = "";
        std::vector<std::string> other_Hs = get_other_H_names(iat, rest);
        if (!other_Hs.empty()) {
            H_at_name_2 = other_Hs[0];
        }

        coot::reduce::add_2_sp3_hydrogens_gemmi(
            H_at_name, H_at_name_2, at_name_1, second_neighb_vec,
            bl_aliph, 109.47, residue_p, true);
        placed_Hs.push_back(H_at_name);
        if (!H_at_name_2.empty()) placed_Hs.push_back(H_at_name_2);

    } else {
        // Default: tetrahedral hydrogen
        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
        std::string at_name_2 = "";
        std::string at_name_3 = "";
        int idx = 0;
        for (unsigned int i2 : iat_neighb_2) {
            if (i2 != iat) {
                if (idx == 0) at_name_2 = rest.atom_info[i2].atom_id_4c;
                else if (idx == 1) at_name_3 = rest.atom_info[i2].atom_id_4c;
                idx++;
            }
        }
        coot::reduce::add_tetrahedral_hydrogen_gemmi(
            H_at_name, at_name_1, at_name_2, at_name_3, "", bl_aliph, *residue_p);
        placed_Hs.push_back(H_at_name);
    }

    return placed_Hs;
}

std::vector<std::string>
place_hydrogen_by_connected_atom_energy_type_gemmi(
    unsigned int iat,
    unsigned int iat_neighb,
    const coot::dictionary_residue_restraints_t& rest,
    gemmi::Residue* residue_p,
    double bl_aliph,
    double bl_arom,
    double bl_amino,
    double bl_oh,
    double bl_sh) {

    std::vector<std::string> v;
    const std::string& energy_type = rest.atom_info[iat_neighb].type_energy;
    const std::string& H_at_name = rest.atom_info[iat].atom_id_4c;
    return place_hydrogen_by_connected_atom_energy_type_gemmi(
        energy_type, iat, iat_neighb, rest, residue_p,
        bl_aliph, bl_arom, bl_amino, bl_oh, bl_sh);
}

} // namespace coot::reduce