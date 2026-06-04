#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bond_length_angle_torsion/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__atoms_with_spinnable_Hs__add__ad3bb6/gemmi/function.hh"

namespace coot {

namespace reduce {

struct torsion_info_t {
    std::string at_name_1;
    std::string at_name_2;
    std::string at_name_3;
    float bond_length;
    float angle_deg;
    float torsion_deg;

    torsion_info_t(const std::string& a1, const std::string& a2, const std::string& a3,
                   float bl, float ang, float tor)
        : at_name_1(a1), at_name_2(a2), at_name_3(a3),
          bond_length(bl), angle_deg(ang), torsion_deg(tor) {}
};

namespace {

static std::string trim_atom_name(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

static const gemmi::Atom* find_atom_in_residue(
    const gemmi::Residue& res,
    const std::string&    name,
    const std::string&    alt_conf)
{
    char alt = alt_conf.empty() ? 0 : alt_conf[0];
    std::string trimmed = trim_atom_name(name);
    for (const auto& a : res.atoms) {
        if (a.name == trimmed && (a.altloc == alt || a.altloc == ' ' || a.altloc == 0)) {
            return &a;
        }
    }
    return nullptr;
}

static gemmi::Atom make_hydrogen_atom(
    const std::string& atom_name,
    const clipper::Coord_orth& pos,
    float bf,
    const std::string& altconf)
{
    gemmi::Atom h;
    h.name = atom_name;
    h.element = gemmi::Element("H");
    h.pos = gemmi::Position(pos.x(), pos.y(), pos.z());
    h.occ = 1.0f;
    h.b_iso = bf;
    if (!altconf.empty() && altconf.length() > 0) {
        h.altloc = altconf[0];
    }
    return h;
}

} // anonymous namespace

inline void add_methyl_Hs_gemmi(
    const std::string &at_name_1,
    const std::string &at_name_2,
    const std::string &at_name_3,
    torsion_info_t torsion_1,
    gemmi::Residue &residue,
    atoms_with_spinnable_Hs &spinables)
{
    std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
    for (unsigned int i = 0; i < alt_confs.size(); i++) {
        clipper::Coord_orth p11;
        clipper::Coord_orth p12;
        clipper::Coord_orth p13;
        bool have_1 = false;

        const gemmi::Atom *at_1 = find_atom_in_residue(residue, torsion_1.at_name_1, alt_confs[i]);
        const gemmi::Atom *at_2 = find_atom_in_residue(residue, torsion_1.at_name_2, alt_confs[i]);
        const gemmi::Atom *at_3 = find_atom_in_residue(residue, torsion_1.at_name_3, alt_confs[i]);

        if (at_1 && at_2 && at_3) {
            have_1 = true;
            p11 = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
                at_1, at_2, at_3,
                torsion_1.bond_length,
                clipper::Util::d2rad(torsion_1.angle_deg),
                clipper::Util::d2rad(torsion_1.torsion_deg));
            p12 = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
                at_1, at_2, at_3,
                torsion_1.bond_length,
                clipper::Util::d2rad(torsion_1.angle_deg),
                clipper::Util::d2rad(torsion_1.torsion_deg + 120));
            p13 = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
                at_1, at_2, at_3,
                torsion_1.bond_length,
                clipper::Util::d2rad(torsion_1.angle_deg),
                clipper::Util::d2rad(torsion_1.torsion_deg - 120));
        }

        if (have_1) {
            clipper::Coord_orth pav_1 = p11;
            clipper::Coord_orth pav_2 = p12;
            clipper::Coord_orth pav_3 = p13;
            float bf = at_3->b_iso;

            gemmi::Atom h1 = make_hydrogen_atom(at_name_1, pav_1, bf, alt_confs[i]);
            gemmi::Atom h2 = make_hydrogen_atom(at_name_2, pav_2, bf, alt_confs[i]);
            gemmi::Atom h3 = make_hydrogen_atom(at_name_3, pav_3, bf, alt_confs[i]);

            residue.atoms.push_back(std::move(h1));
            residue.atoms.push_back(std::move(h2));
            residue.atoms.push_back(std::move(h3));

            gemmi::Atom *new_at_0 = &residue.atoms[residue.atoms.size() - 3];
            gemmi::Atom *new_at_1 = &residue.atoms[residue.atoms.size() - 2];
            gemmi::Atom *new_at_2 = &residue.atoms[residue.atoms.size() - 1];

            std::vector<gemmi::Atom*> h_atoms(3);
            h_atoms[0] = new_at_0;
            h_atoms[1] = new_at_1;
            h_atoms[2] = new_at_2;

            spinables.add_gemmi(
                const_cast<gemmi::Atom*>(at_3),
                atom_with_attached_Hs::hydrogen_t::METHYL,
                h_atoms);
        }
    }
}

} // namespace reduce

} // namespace coot