#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include <utility>
#include <string>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// Reuse connecting_atoms_t_gemmi from the existing port
struct connecting_atoms_t_gemmi {
    connecting_atoms_t_gemmi() : filled_flag(false) {
        upstream_C.first = false;
        upstream_C.second = clipper::Coord_orth(0, 0, 0);
        downstream_N.first = false;
        downstream_N.second = clipper::Coord_orth(0, 0, 0);
    }
    connecting_atoms_t_gemmi(const clipper::Coord_orth &N_pos_in,
                             const clipper::Coord_orth &CA_pos_in,
                             const clipper::Coord_orth &C_pos_in)
        : N_pos(N_pos_in), CA_pos(CA_pos_in), C_pos(C_pos_in), filled_flag(true) {
        upstream_C.first = false;
        upstream_C.second = clipper::Coord_orth(0, 0, 0);
        downstream_N.first = false;
        downstream_N.second = clipper::Coord_orth(0, 0, 0);
    }
    void set_upstream_C(const clipper::Coord_orth &C_pos_in) {
        upstream_C.first = true;
        upstream_C.second = C_pos_in;
    }
    void set_downstream_N(const clipper::Coord_orth &N_pos_in) {
        downstream_N.first = true;
        downstream_N.second = N_pos_in;
    }
    clipper::Coord_orth N_pos;
    clipper::Coord_orth CA_pos;
    clipper::Coord_orth C_pos;
    std::pair<bool, clipper::Coord_orth> upstream_C;
    std::pair<bool, clipper::Coord_orth> downstream_N;
    bool empty() const { return !filled_flag; }
    std::pair<bool, double> get_phi() const {
        if (upstream_C.first) {
            double t = clipper::Coord_orth::torsion(upstream_C.second, N_pos, CA_pos, C_pos);
            return std::pair<bool, double>(true, t);
        } else {
            return std::pair<bool, double>(false, -999.9);
        }
    }
    std::pair<bool, double> get_psi() const {
        if (downstream_N.first) {
            double t = clipper::Coord_orth::torsion(N_pos, CA_pos, C_pos, downstream_N.second);
            return std::pair<bool, double>(true, t);
        } else {
            return std::pair<bool, double>(false, -999.9);
        }
    }
private:
    bool filled_flag;
};

namespace residue_by_phi_psi {

inline connecting_atoms_t_gemmi get_connecting_residue_atoms_gemmi(
    gemmi::Residue* residue_p,
    gemmi::Residue* upstream_neighbour_residue_p,
    gemmi::Residue* downstream_neighbour_residue_p)
{
    connecting_atoms_t_gemmi atoms_in_residue;

    // In gemmi, atom.name is a std::string. For PDB atoms, it's typically
    // a 4-character padded name like " N  ", " C  ", " CA ".
    const auto find_atom = [](gemmi::Residue* res, const std::string& target_name) -> const gemmi::Atom* {
        for (auto& atom : res->atoms) {
            if (atom.name == target_name) return &atom;
        }
        return nullptr;
    };

    const gemmi::Atom* N_at = find_atom(residue_p, "N");
    const gemmi::Atom* C_at = find_atom(residue_p, "C");
    const gemmi::Atom* CA_at = find_atom(residue_p, "CA");

    if (N_at && C_at && CA_at) {
        clipper::Coord_orth N_at_pos = coot::co_gemmi(N_at);
        clipper::Coord_orth CA_at_pos = coot::co_gemmi(CA_at);
        clipper::Coord_orth C_at_pos = coot::co_gemmi(C_at);
        atoms_in_residue = connecting_atoms_t_gemmi(N_at_pos, CA_at_pos, C_at_pos);
    }

    if (upstream_neighbour_residue_p) {
        const gemmi::Atom* up_C = find_atom(upstream_neighbour_residue_p, "C");
        if (up_C) {
            clipper::Coord_orth pos = coot::co_gemmi(up_C);
            atoms_in_residue.set_upstream_C(pos);
        }
    }

    if (downstream_neighbour_residue_p) {
        const gemmi::Atom* down_N = find_atom(downstream_neighbour_residue_p, "N");
        if (down_N) {
            clipper::Coord_orth pos = coot::co_gemmi(down_N);
            atoms_in_residue.set_downstream_N(pos);
        }
    }

    return atoms_in_residue;
}

} // namespace residue_by_phi_psi
} // namespace coot