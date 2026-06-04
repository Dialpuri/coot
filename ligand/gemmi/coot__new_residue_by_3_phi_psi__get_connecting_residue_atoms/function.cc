#include "function.hh"

namespace coot {

connecting_atoms_t_gemmi get_connecting_residue_atoms_gemmi(
    gemmi::Residue* residue_p,
    gemmi::Residue* upstream_neighbour_residue_p,
    gemmi::Residue* downstream_neighbour_residue_p
) {
    connecting_atoms_t_gemmi atoms_in_residue;

    const gemmi::Atom* N_at = nullptr;
    const gemmi::Atom* C_at = nullptr;
    const gemmi::Atom* CA_at = nullptr;

    for (const auto& at : residue_p->atoms) {
        const std::string& atom_name(at.name);
        if (atom_name == " N" || atom_name == "N") N_at = &at;
        if (atom_name == " C" || atom_name == "C") C_at = &at;
        if (atom_name == " CA" || atom_name == "CA") CA_at = &at;
    }

    if (N_at) {
        if (C_at) {
            if (CA_at) {
                clipper::Coord_orth N_at_pos = coot::co_gemmi(N_at);
                clipper::Coord_orth CA_at_pos = coot::co_gemmi(CA_at);
                clipper::Coord_orth C_at_pos = coot::co_gemmi(C_at);
                atoms_in_residue = connecting_atoms_t_gemmi(N_at_pos, CA_at_pos, C_at_pos);
            }
        }
    }

    if (upstream_neighbour_residue_p) {
        for (const auto& at : upstream_neighbour_residue_p->atoms) {
            const std::string& atom_name(at.name);
            if (atom_name == " C" || atom_name == "C") {
                clipper::Coord_orth pos = coot::co_gemmi(&at);
                atoms_in_residue.set_upstream_C(pos);
                break;
            }
        }
    }

    if (downstream_neighbour_residue_p) {
        for (const auto& at : downstream_neighbour_residue_p->atoms) {
            const std::string& atom_name(at.name);
            if (atom_name == " N" || atom_name == "N") {
                clipper::Coord_orth pos = coot::co_gemmi(&at);
                atoms_in_residue.set_downstream_N(pos);
                break;
            }
        }
    }

    return atoms_in_residue;
}

} // namespace coot