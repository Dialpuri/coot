#include "function.hh"
#include <gemmi/calculate.hpp>
#include <cmath>

namespace coot { namespace util {

std::pair<bool, phi_psi_with_residues_t>
get_phi_psi_gemmi(const gemmi::Residue* residue_0,
                  const gemmi::Residue* residue_1,
                  const gemmi::Residue* residue_2) {
    // Check that all residues are valid
    if (!residue_0 || !residue_1 || !residue_2) {
        phi_psi_with_residues_t result;
        return std::make_pair(false, result);
    }

    // Check for consecutive peptide bonding atoms
    const gemmi::Atom* c_prev = nullptr;
    const gemmi::Atom* n_this = nullptr;
    const gemmi::Atom* ca_this = nullptr;
    const gemmi::Atom* c_this = nullptr;
    const gemmi::Atom* n_next = nullptr;

    for (const gemmi::Atom& atom : residue_0->atoms) {
        if (atom.name == "C") c_prev = &atom;
    }
    for (const gemmi::Atom& atom : residue_1->atoms) {
        if (atom.name == "N") n_this = &atom;
        if (atom.name == "CA") ca_this = &atom;
        if (atom.name == "C") c_this = &atom;
    }
    for (const gemmi::Atom& atom : residue_2->atoms) {
        if (atom.name == "N") n_next = &atom;
    }

    // Check that all required atoms are present
    if (!c_prev || !n_this || !ca_this || !c_this || !n_next) {
        phi_psi_with_residues_t result;
        return std::make_pair(false, result);
    }

    // Calculate dihedrals using gemmi::calculate_dihedral
    double phi = gemmi::calculate_dihedral(c_prev->pos, n_this->pos, ca_this->pos, c_this->pos);
    double psi = gemmi::calculate_dihedral(n_this->pos, ca_this->pos, c_this->pos, n_next->pos);

    // Convert from radians to degrees
    phi = phi * 180.0 / M_PI;
    psi = psi * 180.0 / M_PI;

    // Check peptide bonding distances (should be within 2.0A)
    gemmi::Vec3 v1 = c_prev->pos - n_this->pos;
    gemmi::Vec3 v2 = c_this->pos - n_next->pos;
    double dist_1 = v1.length_sq();
    double dist_2 = v2.length_sq();

    bool is_valid = (dist_1 < 4.0) && (dist_2 < 4.0);

    // Create phi_psi_with_residues_t
    phi_psi_with_residues_t result;
    result.phi_ = phi;
    result.psi_ = psi;
    result.residue_prev = residue_0;
    result.residue_this = residue_1;
    result.residue_next = residue_2;

    return std::make_pair(is_valid, result);
}

}} // namespace coot::util