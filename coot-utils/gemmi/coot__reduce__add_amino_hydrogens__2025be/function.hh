#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <iostream>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__is_linked/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bond_length_angle_torsion/gemmi/function.hh"

namespace coot {

namespace reduce {

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
    char alt = alt_conf.empty() ? ' ' : alt_conf[0];
    std::string trimmed = trim_atom_name(name);
    for (const auto& a : res.atoms) {
        if (trim_atom_name(a.name) == trimmed && a.altloc == alt) {
            return &a;
        }
    }
    return nullptr;
}

} // anonymous namespace

inline void add_amino_hydrogens_gemmi(
    const std::string &H_at_name_1,
    const std::string &H_at_name_2,
    const std::string &at_name_1,   // e.g. ND2
    const std::string &at_name_2,   // e.g. CG
    const std::string &at_name_3,   // e.g. OD1
    double bl_amino,
    gemmi::Residue &residue,
    const gemmi::Structure &structure,
    gemmi::Chain *chain)
{
   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
   for (unsigned int i = 0; i < alt_confs.size(); i++) {
      const gemmi::Atom *at_n_1 = find_atom_in_residue(residue, at_name_1, alt_confs[i]);
      const gemmi::Atom *at_n_2 = find_atom_in_residue(residue, at_name_2, alt_confs[i]);
      const gemmi::Atom *at_n_3 = find_atom_in_residue(residue, at_name_3, alt_confs[i]);
      
      if (at_n_1 && at_n_2 && at_n_3) {
         bool add_them = true;
         std::string res_name(residue.name);
         if (res_name == "ASN") {
            gemmi::CRA cra{chain, &residue, const_cast<gemmi::Atom*>(at_n_1)};
            if (coot::reduce::is_linked_gemmi(at_name_1, cra, structure)) {
               // Don't add hydrogens on ND2 linked to NAG
               add_them = false;
            }
         }
         if (add_them) {
            clipper::Coord_orth Hp1 = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
                at_n_3, at_n_2, at_n_1,
                bl_amino,
                clipper::Util::d2rad(120),
                clipper::Util::d2rad(180));
            clipper::Coord_orth Hp2 = coot::reduce::position_by_bond_length_angle_torsion_gemmi(
                at_n_3, at_n_2, at_n_1,
                bl_amino,
                clipper::Util::d2rad(120),
                clipper::Util::d2rad(0));
            float bf = static_cast<float>(at_n_1->b_iso);
            coot::reduce::add_hydrogen_atom_gemmi(H_at_name_1, Hp1, bf, alt_confs[i], residue);
            coot::reduce::add_hydrogen_atom_gemmi(H_at_name_2, Hp2, bf, alt_confs[i], residue);
         }
      } else {
         std::cout << "Fail Residue " << residue.name
                   << " alt-conf \"" << alt_confs[i] << "\""
                   << " failed in add_amino_hydrogens" << std::endl;
      }
   }
}

} // namespace reduce

} // namespace coot