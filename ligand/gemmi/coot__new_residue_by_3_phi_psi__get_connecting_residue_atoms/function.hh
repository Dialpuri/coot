#pragma once

#include <clipper/core/xmap.h>
#include <gemmi/model.hpp>
#include <utility>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

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

connecting_atoms_t_gemmi get_connecting_residue_atoms_gemmi(
    gemmi::Residue* residue_p,
    gemmi::Residue* upstream_neighbour_residue_p,
    gemmi::Residue* downstream_neighbour_residue_p
);

} // namespace coot