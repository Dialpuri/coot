#include "function.hh"
#include <algorithm>

namespace coot {
namespace rotamer_helper {

static void swap_atom_names(
    gemmi::Residue* res,
    const std::vector<std::pair<std::string, std::string>>& swapper_atoms) {
    for (unsigned int iswap = 0; iswap < swapper_atoms.size(); iswap++) {
        bool ifound = false;
        for (size_t ifirst = 0; ifirst < res->atoms.size() && !ifound; ifirst++) {
            if (res->atoms[ifirst].name == swapper_atoms[iswap].first) {
                for (size_t isec = 0; isec < res->atoms.size(); isec++) {
                    if (res->atoms[isec].name == swapper_atoms[iswap].second) {
                        std::swap(res->atoms[ifirst].name, res->atoms[isec].name);
                        ifound = true;
                    }
                    if (ifound) break;
                }
            }
        }
    }
}

static bool has_any_swap_pair(
    const gemmi::Residue* res,
    const std::vector<std::pair<std::string, std::string>>& swapper_atoms) {
    for (unsigned int iswap = 0; iswap < swapper_atoms.size(); iswap++) {
        bool first_found = false;
        bool sec_found = false;
        for (const auto& atom : res->atoms) {
            if (atom.name == swapper_atoms[iswap].first) first_found = true;
            if (atom.name == swapper_atoms[iswap].second) sec_found = true;
        }
        if (first_found && sec_found) return true;
    }
    return false;
}

int optimize_rotamer_by_atom_names_gemmi(gemmi::CRA cra, bool apply_swap_if_found) {
    int imoved = 0;

    if (!cra.residue)
        return imoved;

    // Get initial rotamer probability
    coot::rotamer_probability_info_t p_init =
        coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);

    std::string residue_name(cra.residue->name);

    if ((residue_name == "PHE") ||
        (residue_name == "ASP") ||
        (residue_name == "GLU") ||
        (residue_name == "TYR")) {

        std::vector<std::pair<std::string, std::string>> swapper_atoms;

        if (residue_name == "PHE") {
            swapper_atoms.push_back(std::make_pair("CD1", "CD2"));
            swapper_atoms.push_back(std::make_pair("CE1", "CE2"));
        }
        if (residue_name == "ASP") {
            swapper_atoms.push_back(std::make_pair("OD1", "OD2"));
        }
        if (residue_name == "GLU") {
            swapper_atoms.push_back(std::make_pair("OE1", "OE2"));
        }
        if (residue_name == "TYR") {
            swapper_atoms.push_back(std::make_pair("CD1", "CD2"));
            swapper_atoms.push_back(std::make_pair("CE1", "CE2"));
        }

        // Check if this residue has the swap pairs
        if (!has_any_swap_pair(cra.residue, swapper_atoms))
            return imoved;

        // If initial probability state is bad, we can't compare
        if (!p_init.state)
            return imoved;

        // Swap atom names in-place (this is the "try swap" step)
        swap_atom_names(cra.residue, swapper_atoms);

        // Get probability of the swapped configuration
        coot::rotamer_probability_info_t p_swapped =
            coot::rotamer_helper::probability_of_this_rotamer_gemmi(cra);

        // Compare probabilities
        if (p_swapped.state) {
            if (p_swapped.probability > p_init.probability) {
                // Swapped is better
                if (apply_swap_if_found) {
                    imoved = 1;
                    // Keep swapped atoms in place
                } else {
                    // Swap back to original
                    swap_atom_names(cra.residue, swapper_atoms);
                }
            } else {
                // Original was better; swap back
                swap_atom_names(cra.residue, swapper_atoms);
            }
        } else {
            // Swapped probability is bad; swap back to original
            swap_atom_names(cra.residue, swapper_atoms);
        }
    }

    return imoved;
}

} // namespace rotamer_helper
} // namespace coot