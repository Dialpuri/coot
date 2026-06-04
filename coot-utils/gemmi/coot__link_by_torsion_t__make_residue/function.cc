#include "function.hh"
#include <gemmi/elem.hpp>
#include "coot/coot-utils/coot-coord-utils.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_by_torsion_t__pos/gemmi/function.hh"

namespace coot {

gemmi::Residue* make_residue_gemmi(
    const std::string& new_residue_type,
    int new_res_no,
    const std::vector<atom_torsion_data_t>& atom_torsions,
    const gemmi::Residue& base_residue,
    float b_factor) {

   gemmi::Residue* r = nullptr;
   if (!atom_torsions.empty()) {
      r = new gemmi::Residue;
      r->name = new_residue_type;
      r->seqid.num.value = new_res_no;

      for (unsigned int i = 0; i < atom_torsions.size(); i++) {
         const atom_torsion_data_t& gat = atom_torsions[i];

         // filled_atom_name — inline (no _gemmi port)
         std::string f = "XXXX";
         std::string::size_type len = gat.atom_name.length();
         if (len == 1) f = " " + gat.atom_name + "  ";
         if (len == 2) f = " " + gat.atom_name + " ";
         if (len == 3) f = " " + gat.atom_name;

         // Use the verified _gemmi port for pos
         clipper::Coord_orth p = coot::pos_gemmi(
             gat.atom_name,
             gat.prior_atom_1,
             gat.prior_atom_2,
             gat.prior_atom_3,
             base_residue,
             *r,
             gat.bond_length,
             gat.angle_deg,
             gat.torsion_deg);

         gemmi::Atom atom;
         atom.name = f;
         atom.element = gemmi::Element(gat.element);
         atom.pos = gemmi::Position(p.x(), p.y(), p.z());
         atom.occ = 1.0;
         atom.b_iso = b_factor;
         r->atoms.push_back(atom);
      }
   }
   return r;
}

} // namespace coot