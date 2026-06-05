#pragma once

#include <gemmi/model.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// Helper: pad atom name to 4 characters matching MMDB GetAtomName() convention
// PDB convention: 1-char names get leading space (" N  "), others are right-padded ("CB  ")
static std::string pad_atom_name(const std::string &name) {
  std::string s = name;
  if (s.size() == 1) {
    s = " " + s;  // single-char atoms: " N  "
  }
  while (s.size() < 4) s += ' ';
  return s;
}

// Helper: format atom_spec_t equivalent for a gemmi atom
// Matches original: "[spec: model N \"C\"    R \"I\" \" AN \" \"A\"]"
inline std::string atom_spec_string(
  const gemmi::Atom &at,
  int model_number,
  const std::string &chain_id,
  int res_no,
  const std::string &ins_code)
{
  std::ostringstream s;
  std::string alt_loc;
  if (at.altloc != '\0' && at.altloc != ' ') {
    alt_loc = std::string(1, at.altloc);
  }
  s << "[spec: model " << model_number << " \"" << chain_id << "\" "
    << std::setw(4) << std::right << res_no << " \""
    << ins_code << "\" \""
    << pad_atom_name(at.name) << "\" \""
    << alt_loc << "\"]";
  return s.str();
}

// gemmi port of coot::glyco_tree_t::output_internal_distances
// Outputs atom-atom distance information between a central residue and
// a list of residues, within a distance criterion.
inline void output_internal_distances_gemmi(
  const gemmi::Residue &residue,
  const std::vector<const gemmi::Residue *> residues,
  double dist_crit,
  std::ofstream &f,
  int model_number,
  const std::string &chain_id)
{
  double dist_min = 2.66;

  // Collect atoms from the central residue
  std::vector<const gemmi::Atom *> central_atoms;
  for (const auto &at : residue.atoms) {
    central_atoms.push_back(&at);
  }

  // Get residue sequence number and insertion code from the central residue
  int res_no = residue.seqid.num.value;
  std::string ins_code;
  char icode = residue.seqid.icode;
  if (icode != ' ') ins_code = std::string(1, icode);

  for (const auto *res_p : residues) {
    if (res_p != &residue) {
      int res_no_at = res_p->seqid.num.value;
      std::string ins_code_at;
      char icode_at = res_p->seqid.icode;
      if (icode_at != ' ') ins_code_at = std::string(1, icode_at);

      for (const auto &at : res_p->atoms) {
        // Skip atoms with all-zero coordinates (terminated atoms)
        if (at.pos.x == 0 && at.pos.y == 0 && at.pos.z == 0) continue;

        clipper::Coord_orth pos_atom = co_gemmi(&at);

        for (const auto *at_c : central_atoms) {
          if (at_c != &at) {
            // Skip terminated atoms
            if (at_c->pos.x == 0 && at_c->pos.y == 0 && at_c->pos.z == 0) continue;

            clipper::Coord_orth pos_central_atom = co_gemmi(at_c);
            double d = clipper::Coord_orth::length(pos_atom, pos_central_atom);

            if (d < dist_crit && d > dist_min) {
              std::string spec_at = atom_spec_string(
                at, model_number, chain_id, res_no_at, ins_code_at);
              std::string spec_c = atom_spec_string(
                *at_c, model_number, chain_id, res_no, ins_code);

              f << " " << spec_at << " " << spec_c << " " << d << std::endl;
            }
          }
        }
      }
    }
  }
}

} // namespace coot