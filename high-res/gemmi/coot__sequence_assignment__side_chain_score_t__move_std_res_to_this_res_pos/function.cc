#include "function.hh"

void
coot::sequence_assignment::side_chain_score_t::move_std_res_to_this_res_pos_gemmi(const clipper::RTop_orth &rtop,
                                                                                  gemmi::Residue *std_residue) {
  for (gemmi::Atom &atom : std_residue->atoms) {
    clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
    clipper::Coord_orth rotted = co.transform(rtop);
    atom.pos.x = rotted.x();
    atom.pos.y = rotted.y();
    atom.pos.z = rotted.z();
  }
}