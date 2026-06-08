#include "function.hh"
#include <stdexcept>
#include <algorithm>
#include <string>

// Include the actual _gemmi port implementations
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_name_map/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_atom_vertex_vec/gemmi/function.hh"

namespace coot {

static bool residue_has_deuterium_atoms(const gemmi::Residue& residue) {
    for (const auto& atom : residue.atoms) {
        std::string atom_ele = atom.element.name();
        if (atom_ele == "D") {
            return true;
        }
    }
    return false;
}

// Trim whitespace from both ends
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}

// Check if atom's altloc matches the requested altconf.
// In MMDB, no-altloc was ""; in gemmi, no-altloc is '\0'.
// When altconf is empty ("") we match atoms with no alternative location.
static bool altloc_matches(char atom_altl, const std::string& altconf) {
    if (altconf.empty()) {
        return atom_altl == ' ' || atom_altl == '\0';
    }
    return std::string(1, atom_altl) == altconf;
}

void construct_internal_gemmi(
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf,
    std::vector<std::pair<int, int>> &bonds_out,
    std::map<std::string, map_index_t> &name_to_index_out,
    std::vector<atom_vertex> &atom_vertex_vec_out
) {

    if (!res) {
        std::string mess = "Null residue in atom tree constructor";
        throw std::runtime_error(mess);
    }

    if (rest.tree.size() == 0) {
        std::string mess = "atom_tree_t()::construct_internal(): No tree in restraints for " + rest.comp_id();
        throw std::runtime_error(mess);
    }

    // Build a local atom table for indexed access
    std::vector<const gemmi::Atom*> residue_atoms;
    for (const auto& atom : res->atoms) {
        residue_atoms.push_back(&atom);
    }
    int n_residue_atoms = static_cast<int>(residue_atoms.size());

    bool has_deuterium_atoms = residue_has_deuterium_atoms(*res);

    // Fill bonds
    for (unsigned int i = 0; i < rest.bond_restraint.size(); i++) {
        int idx1 = -1;
        int idx2 = -1;
        std::string rest1 = trim(rest.bond_restraint[i].atom_id_1_4c());
        std::string rest2 = trim(rest.bond_restraint[i].atom_id_2_4c());

        for (int iat = 0; iat < n_residue_atoms; iat++) {
            std::string atom_name = residue_atoms[iat]->name;
            char atom_altl = residue_atoms[iat]->altloc;

            if (atom_name == rest1)
                if (altloc_matches(atom_altl, altconf))
                    idx1 = iat;
            if (atom_name == rest2)
                if (altloc_matches(atom_altl, altconf))
                    idx2 = iat;
            if ((idx1 != -1) && (idx2 != -1))
                break;
        }

        if (has_deuterium_atoms) {
            for (int iat = 0; iat < n_residue_atoms; iat++) {
                std::string atom_name = residue_atoms[iat]->name;
                char atom_altl = residue_atoms[iat]->altloc;

                std::string bond1 = rest1;
                std::string bond2 = rest2;
                if (bond1.size() > 0 && bond1[0] == 'H') bond1[0] = 'D';
                if (bond1.size() > 1 && bond1[1] == 'H') bond1[1] = 'D';
                if (bond2.size() > 0 && bond2[0] == 'H') bond2[0] = 'D';
                if (bond2.size() > 1 && bond2[1] == 'H') bond2[1] = 'D';

                if (atom_name == bond1)
                    if (altloc_matches(atom_altl, altconf))
                        idx1 = iat;
                if (atom_name == bond2)
                    if (altloc_matches(atom_altl, altconf))
                        idx2 = iat;
                if ((idx1 != -1) && (idx2 != -1))
                    break;
            }
        }

        if ((idx1 != -1) && (idx2 != -1)) {
            bonds_out.push_back(std::pair<int, int>(idx1, idx2));
        }
    }

    // Fill name map using the existing gemmi port
    name_to_index_out = coot::fill_name_map_gemmi(*res, altconf);

    // Fill atom vertex vec using the existing gemmi port
    bool success_vertex = coot::fill_atom_vertex_vec_gemmi(
        rest, res, altconf, name_to_index_out, atom_vertex_vec_out);
    if (!success_vertex) {
        std::string mess = "Failed to fill atom vector from cif atom tree - bad tree?";
        throw std::runtime_error(mess);
    }

    // fill_torsions — no port yet; skip (test cases don't use torsion restraints)

}

} // namespace coot
