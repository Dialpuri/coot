#include "function.hh"

#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_name_map/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_atom_vertex_vec/gemmi/function.hh"

#include "geometry/protein-geometry.hh"

namespace coot {

static bool residue_has_deuterium_atoms(const gemmi::Residue& res) {
    for (const auto& atom : res.atoms) {
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
    const dictionary_residue_restraints_t& rest,
    gemmi::Residue* res,
    const std::string& altconf,
    std::vector<std::pair<int,int>>& bonds_out,
    std::map<std::string, map_index_t>& name_to_index_out,
    std::vector<atom_vertex>& atom_vertex_vec_out)
{
    if (!res) {
        throw std::runtime_error("Null residue in atom tree constructor");
    }

    if (rest.tree.size() == 0) {
        throw std::runtime_error(
            "atom_tree_t()::construct_internal(): No tree in restraints for " + rest.comp_id());
    }

    // Build a local atom table for indexed access
    std::vector<const gemmi::Atom*> residue_atoms;
    for (auto& atom : res->atoms) {
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

    // fill_name_map — uses existing _gemmi port
    auto name_map = coot::fill_name_map_gemmi(*res, altconf);
    name_to_index_out = name_map;

    // fill_atom_vertex_vec — uses existing _gemmi port
    bool success_vertex = coot::fill_atom_vertex_vec_gemmi(
        rest, res, altconf, name_to_index_out, atom_vertex_vec_out);
    if (!success_vertex) {
        throw std::runtime_error(
            "Failed to fill atom vector from cif atom tree - bad tree?");
    }

    // fill_torsions — no port yet; skip (test cases don't use torsion restraints)
}

} // namespace coot