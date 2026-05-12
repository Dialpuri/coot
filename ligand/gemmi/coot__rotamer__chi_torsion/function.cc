#include "function.hh"

namespace coot {
namespace rotamer {

double chi_torsion_gemmi(const std::vector<int> &chi_angle_atoms_indices,
                         const std::vector<gemmi::Atom>& residue_atoms) {
   double tors = 0.0;

   std::vector<clipper::Coord_orth> a;
   
   for (unsigned int ich_at=0; ich_at<chi_angle_atoms_indices.size(); ich_at++) {
      const gemmi::Atom& atom = residue_atoms[chi_angle_atoms_indices[ich_at]];
      a.push_back(clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
   }

   double ctorsion = clipper::Coord_orth::torsion(a[0], a[1], a[2], a[3]);
   tors = clipper::Util::rad2d(ctorsion);

   return tors;
}

} // namespace rotamer
} // namespace coot