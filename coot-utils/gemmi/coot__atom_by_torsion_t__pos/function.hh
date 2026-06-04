#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <string>
#include <utility>
#include <stdexcept>
#include <sstream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

static inline std::string int_to_string(int n) {
    std::ostringstream os;
    os << n;
    return os.str();
}

/// Gemmi port of coot::atom_by_torsion_t::pos
///
/// Computes a new atom position from three reference atoms, a bond length,
/// an angle (degrees) and a torsion (degrees).
///
/// prior_atom_X.first == true  → look in base_res
/// prior_atom_X.first == false → look in ext_res
///
clipper::Coord_orth pos_gemmi(
    const std::string& atom_name,
    std::pair<bool, std::string> prior_atom_1,
    std::pair<bool, std::string> prior_atom_2,
    std::pair<bool, std::string> prior_atom_3,
    const gemmi::Residue& base_res,
    const gemmi::Residue& ext_res,
    double bond_length,
    double angle_deg,
    double torsion_deg)
{
    const gemmi::Atom* at_1 = nullptr;
    const gemmi::Atom* at_2 = nullptr;
    const gemmi::Atom* at_3 = nullptr;

    auto get_atom = [](const gemmi::Residue& from_base,
                       const gemmi::Residue& from_ext,
                       bool use_base,
                       const std::string& atom_name) -> const gemmi::Atom* {
        const gemmi::Residue& res = use_base ? from_base : from_ext;
        return res.find_atom(atom_name, '*');
    };

    at_1 = get_atom(base_res, ext_res, prior_atom_1.first, prior_atom_1.second);
    at_2 = get_atom(base_res, ext_res, prior_atom_2.first, prior_atom_2.second);
    at_3 = get_atom(base_res, ext_res, prior_atom_3.first, prior_atom_3.second);

    if (at_1 && at_2 && at_3) {
        clipper::Coord_orth p1 = coot::co_gemmi(at_1);
        clipper::Coord_orth p2 = coot::co_gemmi(at_2);
        clipper::Coord_orth p3 = coot::co_gemmi(at_3);
        clipper::Coord_orth new_pos = clipper::Coord_orth(
            p3, p2, p1,
            bond_length,
            clipper::Util::d2rad(angle_deg),
            clipper::Util::d2rad(torsion_deg));
        return new_pos;
    } else {
        unsigned int n_missing = 0;
        if (!at_1) n_missing++;
        if (!at_2) n_missing++;
        if (!at_3) n_missing++;

        std::string m = "missing atom";
        if (n_missing > 1)
            m += "s";
        m += " in atom_by_torsion_t::pos() when positioning ";
        m += atom_name;
        m += " : ";
        if (!at_1) m += " at_1: " + prior_atom_1.second;
        if (!at_2) m += " at_2: " + prior_atom_2.second;
        if (!at_3) m += " at_3: " + prior_atom_3.second;
        m += " of ";
        m += int_to_string(static_cast<int>(base_res.atoms.size()));
        m += " base atoms";
        throw std::runtime_error(m);
    }
}

} // namespace coot