#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <string>
#include <vector>
#include <iostream>

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

namespace coot {
namespace nomenclature {

int test_and_fix_PHE_TYR_nomenclature_errors_gemmi(gemmi::CRA cra, bool apply_swap_if_found) {
  if (!cra.residue) return 0;

  int iswapped = 0;

  // Collect unique alt confs from CD1 atoms
  std::vector<char> alt_conf_list;
  for (auto& atom : cra.residue->atoms) {
    if (atom.name == "CD1") {
      bool already = false;
      for (char ac : alt_conf_list) {
        if (ac == atom.altloc) { already = true; break; }
      }
      if (!already) alt_conf_list.push_back(atom.altloc);
    }
  }

  if (alt_conf_list.empty()) return 0;

  for (unsigned int ialtconf = 0; ialtconf < alt_conf_list.size(); ialtconf++) {
    gemmi::Atom* CA  = nullptr;
    gemmi::Atom* CB  = nullptr;
    gemmi::Atom* CG  = nullptr;
    gemmi::Atom* CD1 = nullptr;
    gemmi::Atom* CD2 = nullptr;

    char alt = alt_conf_list[ialtconf];

    // First pass: find atoms with matching alt conf
    for (auto& atom : cra.residue->atoms) {
      if (atom.altloc == alt) {
        if (atom.name == "CA") CA  = &atom;
        if (atom.name == "CB") CB  = &atom;
        if (atom.name == "CG") CG  = &atom;
        if (atom.name == "CD1") CD1 = &atom;
        if (atom.name == "CD2") CD2 = &atom;
      }
    }

    // Fallback: if CA, CB, or CG missing, try alt conf '\0' (default)
    if (!CA || !CB || !CG) {
      for (auto& atom : cra.residue->atoms) {
        if (atom.altloc == '\0') {
          if (!CA && atom.name == "CA") CA  = &atom;
          if (!CB && atom.name == "CB") CB  = &atom;
          if (!CG && atom.name == "CG") CG  = &atom;
          if (!CD2 && atom.name == "CD2") CD2 = &atom;
        }
      }
    }

    if (CA && CB && CG && CD1) {
      // gemmi::calculate_dihedral returns radians; convert to degrees
      double tors = gemmi::calculate_dihedral(CA->pos, CB->pos, CG->pos, CD1->pos);
      tors = tors * 180.0 / M_PI;

      if (tors < -90.0 || tors > 90.0) {
        if (CD2) {
          double cg2_tors = gemmi::calculate_dihedral(CA->pos, CB->pos, CG->pos, CD2->pos);
          cg2_tors = cg2_tors * 180.0 / M_PI;

          if (cg2_tors > -90.0 && cg2_tors < 90.0) {
            gemmi::Atom* CE1 = nullptr;
            gemmi::Atom* CE2 = nullptr;
            gemmi::Atom* HD1 = nullptr;
            gemmi::Atom* HD2 = nullptr;
            gemmi::Atom* HE1 = nullptr;
            gemmi::Atom* HE2 = nullptr;
            for (auto& atom : cra.residue->atoms) {
              if (atom.altloc == alt) {
                if (atom.name == "CE1") CE1 = &atom;
                if (atom.name == "CE2") CE2 = &atom;
                if (atom.name == "HD1") HD1 = &atom;
                if (atom.name == "HD2") HD2 = &atom;
                if (atom.name == "HE1") HE1 = &atom;
                if (atom.name == "HE2") HE2 = &atom;
              }
            }

            if (CE1 && CE2) {
              if (apply_swap_if_found) {
                gemmi::Position pos_cd1 = CD1->pos;
                gemmi::Position pos_cd2 = CD2->pos;
                gemmi::Position pos_ce1 = CE1->pos;
                gemmi::Position pos_ce2 = CE2->pos;
                CD1->pos = pos_cd2;
                CD2->pos = pos_cd1;
                CE1->pos = pos_ce2;
                CE2->pos = pos_ce1;
              }
              iswapped++;
            }

            if (HD1 && HD2 && apply_swap_if_found) {
              HD1->name = "HD2";
              HD2->name = "HD1";
            }
            if (HE1 && HE2 && apply_swap_if_found) {
              HE1->name = "HE2";
              HE2->name = "HE1";
            }
          }
        }
      }
    } else {
      std::cout << "WARNING:: PHE/TYR nomenc: missing atoms in "
                << (cra.chain ? std::string(1, cra.chain->name[0]) : "?") << " "
                << cra.residue->seqid.num.value << " " << cra.residue->name
                << std::endl;
    }
  }

  return iswapped;
}

} // namespace nomenclature
} // namespace coot