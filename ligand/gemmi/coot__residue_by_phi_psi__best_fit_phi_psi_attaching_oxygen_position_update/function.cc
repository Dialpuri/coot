#include "function.hh"
#include "coot/mini-mol/mini-mol-utils.hh"

clipper::Coord_orth coot::residue_by_phi_psi::best_fit_phi_psi_attaching_oxygen_position_update_gemmi(
    const minimol::molecule &mm, const gemmi::Residue *residue_p) {

  clipper::Coord_orth pos(0, 0, 0);

  // Build a minimol::residue from the gemmi residue
  minimol::residue res_with_CA_C;
  res_with_CA_C.seqnum = residue_p->seqid.num.value;
  res_with_CA_C.ins_code = std::string(1, residue_p->seqid.icode);
  res_with_CA_C.name = residue_p->name;
  for (const gemmi::Atom &gemmi_atom : residue_p->atoms) {
    res_with_CA_C.addatom(gemmi_atom.name,
                          gemmi_atom.element.name(),
                          clipper::Coord_orth(gemmi_atom.pos.x,
                                              gemmi_atom.pos.y,
                                              gemmi_atom.pos.z),
                          "",
                          gemmi_atom.b_iso,
                          gemmi_atom.occ);
  }

  if (mm.fragments.size() == 1) {
    const minimol::fragment &frag = mm.fragments[0];
    int idx_first = frag.first_residue();
    const minimol::residue &res_with_N = frag[idx_first];

    // Call o_position (free function in coot namespace from mini-mol-utils.hh)
    std::pair<bool, clipper::Coord_orth> o_pos =
        coot::o_position(res_with_CA_C, res_with_N);
    if (o_pos.first) {
      pos = o_pos.second;
    }
  }

  return pos;
}