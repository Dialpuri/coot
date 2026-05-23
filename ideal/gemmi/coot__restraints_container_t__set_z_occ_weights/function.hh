#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>
#include <iostream>

namespace coot {

struct restraints_container_t {
    std::vector<gemmi::Atom> atom;
    size_t n_atoms;
    std::vector<float> atom_z_occ_weight;
    bool do_neutron_refinement;
    bool cryo_em_mode;
    
    void set_z_occ_weights_gemmi();
};

// Helper functions from coot::util namespace (no _gemmi ports available)
inline std::vector<std::pair<std::string, int>> atomic_number_atom_list() {
    return {
        {"H", 1}, {"D", 1}, {"T", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
        {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
        {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
        {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
        {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42},
        {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
        {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54},
        {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
        {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66},
        {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
        {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78},
        {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84},
        {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
        {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96},
        {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102},
        {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107},
        {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112},
        {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}
    };
}

inline int atomic_number(const std::string& atom_name, const std::vector<std::pair<std::string, int>>& atom_list) {
    std::string el = atom_name;
    if (el.size() >= 2) {
        if (el[0] >= 'a' && el[0] <= 'z') el[0] = el[0] - 'a' + 'A';
        if (el.size() > 1 && el[1] >= 'A' && el[1] <= 'Z') el[1] = el[1] - 'A' + 'a';
    } else if (el.size() == 1 && el[0] >= 'a' && el[0] <= 'z') {
        el[0] = el[0] - 'a' + 'A';
    }
    
    for (const auto& p : atom_list) {
        if (p.first == el) {
            return p.second;
        }
    }
    return -1;
}

inline double neutron_occupancy(const std::string& element, int formal_charge) {
    if (element == "H" || element == "D" || element == "T") return 1.0;
    if (element == "N") return 1.0;
    if (element == "O") return 1.0;
    if (element == "C") return 1.0;
    if (element == "S") return 1.0;
    if (element == "P") return 1.0;
    return 1.0;
}

inline bool is_main_chain_or_cb_p(const gemmi::Atom& atom) {
    const std::string& name = atom.name;
    if (name == "CA" || name == "C" || name == "N" || name == "O") return true;
    if (name == "CB") return true;
    return false;
}

inline void restraints_container_t::set_z_occ_weights_gemmi() {
    atom_z_occ_weight.resize(n_atoms);
    std::vector<std::pair<std::string, int>> atom_list = atomic_number_atom_list();
    for (size_t i = 0; i < n_atoms; i++) {
        const gemmi::Atom& at = atom[i];
        if (at.element != gemmi::Element("")) {
            std::string element = at.element.name();
            double z = atomic_number(element, atom_list);
            double weight = 1.0;
            double occupancy = at.occ;
            if (occupancy > 1.0) occupancy = 1.0;
            if (do_neutron_refinement) {
                int formal_charge = 0;
                occupancy = neutron_occupancy(element, formal_charge);
            }
            if (cryo_em_mode) {
                if (!is_main_chain_or_cb_p(at)) {
                    weight = 0.2;
                }
                std::string at_name = at.name;
                if (at_name == " O  ") {
                    weight = 0.4;
                }
            }

            if (z < 0.0) {
                std::cout << "WARNING:: init_shared_post() atom " << i << " " 
                          << "Unknown element \"" << element << "\"" << std::endl;
                z = 6.0;
            }
            atom_z_occ_weight[i] = static_cast<float>(weight * z * occupancy);
        }
    }
}

} // namespace coot