#include "function.hh"

void coot::util::transform_selection_gemmi(
    gemmi::Model& model,
    const std::string& chain_name,
    const clipper::RTop_orth& rtop)
{
  for (gemmi::Chain& chain : model.chains) {
    if (chain.name != chain_name)
      continue;
    for (gemmi::Residue& res : chain.residues) {
      for (gemmi::Atom& atom : res.atoms) {
        clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
        clipper::Coord_orth trans_pos = co.transform(rtop);
        atom.pos.x = trans_pos.x();
        atom.pos.y = trans_pos.y();
        atom.pos.z = trans_pos.z();
      }
    }
  }
}