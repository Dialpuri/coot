#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/math.hpp>
#include <memory>
#include <string>
#include <vector>
#include "mini-mol/mini-mol.hh"

namespace coot {
namespace minimol {

/// gemmi port of coot::minimol::residue::make_residue()
/// Converts a mini-mol residue into a gemmi::Residue.
/// Returns nullptr if the residue has no atoms.
inline std::unique_ptr<gemmi::Residue> make_residue_gemmi(coot::minimol::residue const& res) {
  if (res.atoms.empty()) return nullptr;

  auto gemmi_res = std::make_unique<gemmi::Residue>();

  // Set residue identity
  gemmi_res->name = res.name;
  gemmi_res->seqid.num.value = res.seqnum;
  // Insertion code: gemmi uses ' ' for "no insertion code", mini-mol uses ""
  if (res.ins_code.empty()) {
    gemmi_res->seqid.icode = ' ';
  } else {
    gemmi_res->seqid.icode = res.ins_code[0];
  }

  // Populate atoms
  for (unsigned int iat = 0; iat < res.atoms.size(); iat++) {
    gemmi::Atom atom;
    atom.name = res.atoms[iat].name;
    atom.element = gemmi::Element(res.atoms[iat].element);
    atom.pos = gemmi::Position(res.atoms[iat].pos.x(),
                               res.atoms[iat].pos.y(),
                               res.atoms[iat].pos.z());
    atom.occ = static_cast<double>(res.atoms[iat].occupancy);
    atom.b_iso = static_cast<double>(res.atoms[iat].temperature_factor);
    if (!res.atoms[iat].altLoc.empty()) {
      atom.altloc = res.atoms[iat].altLoc[0];
    } else {
      atom.altloc = '\0';
    }
    gemmi_res->atoms.push_back(atom);
  }

  return gemmi_res;
}

} // namespace minimol
} // namespace coot