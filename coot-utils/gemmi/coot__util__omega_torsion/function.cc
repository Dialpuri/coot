#include "function.hh"

namespace coot { namespace util {

std::pair<bool, double>
omega_torsion_gemmi(const gemmi::Residue* C_residue,
                    const gemmi::Residue* N_residue,
                    const std::string& altconf) {
    if (!C_residue) return std::make_pair(false, 0.0);
    if (!N_residue) return std::make_pair(false, 0.0);

    const gemmi::Atom* C_residue_CA_atom = nullptr;
    const gemmi::Atom* C_residue_C_atom = nullptr;
    
    for (const auto& atom : C_residue->atoms) {
        if (atom.name == "CA") {
            if (altconf.empty()) {
                if (atom.altloc == '\0' || atom.altloc == ' ') {
                    C_residue_CA_atom = &atom;
                }
            } else {
                if (atom.altloc == altconf[0]) {
                    C_residue_CA_atom = &atom;
                }
            }
        }
        if (atom.name == "C") {
            if (altconf.empty()) {
                if (atom.altloc == '\0' || atom.altloc == ' ') {
                    C_residue_C_atom = &atom;
                }
            } else {
                if (atom.altloc == altconf[0]) {
                    C_residue_C_atom = &atom;
                }
            }
        }
    }

    const gemmi::Atom* N_residue_CA_atom = nullptr;
    const gemmi::Atom* N_residue_N_atom = nullptr;
    
    for (const auto& atom : N_residue->atoms) {
        if (atom.name == "CA") {
            if (altconf.empty()) {
                if (atom.altloc == '\0' || atom.altloc == ' ') {
                    N_residue_CA_atom = &atom;
                }
            } else {
                if (atom.altloc == altconf[0]) {
                    N_residue_CA_atom = &atom;
                }
            }
        }
        if (atom.name == "N") {
            if (altconf.empty()) {
                if (atom.altloc == '\0' || atom.altloc == ' ') {
                    N_residue_N_atom = &atom;
                }
            } else {
                if (atom.altloc == altconf[0]) {
                    N_residue_N_atom = &atom;
                }
            }
        }
    }

    if (C_residue_CA_atom && C_residue_C_atom && N_residue_N_atom && N_residue_CA_atom) {
        clipper::Coord_orth ca1(C_residue_CA_atom->pos.x,
                                C_residue_CA_atom->pos.y,
                                C_residue_CA_atom->pos.z);
        clipper::Coord_orth c1(C_residue_C_atom->pos.x,
                               C_residue_C_atom->pos.y,
                               C_residue_C_atom->pos.z);
        clipper::Coord_orth ca2(N_residue_CA_atom->pos.x,
                                N_residue_CA_atom->pos.y,
                                N_residue_CA_atom->pos.z);
        clipper::Coord_orth n2(N_residue_N_atom->pos.x,
                               N_residue_N_atom->pos.y,
                               N_residue_N_atom->pos.z);

        double omega = clipper::Coord_orth::torsion(ca1, c1, n2, ca2);
        return std::make_pair(true, omega);
    }
    
    return std::make_pair(false, 0.0);
}

}} // namespace coot::util