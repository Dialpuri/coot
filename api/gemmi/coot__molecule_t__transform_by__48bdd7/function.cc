#include "function.hh"

namespace coot {

void transform_by_gemmi(gemmi::Structure& st, std::array<std::array<double, 4>, 4> mat) {
   // Extract 3x3 rotation matrix and translation vector from 4x4 matrix
   clipper::Mat33<double> clipper_mat(
       mat[0][0], mat[0][1], mat[0][2],
       mat[1][0], mat[1][1], mat[1][2],
       mat[2][0], mat[2][1], mat[2][2]
   );
   clipper::Coord_orth trans_vec(mat[0][3], mat[1][3], mat[2][3]);
   clipper::RTop_orth rtop(clipper_mat, trans_vec);

   // Iterate over all atoms in the structure and apply transformation
   for (gemmi::Model &model : st.models) {
      for (gemmi::Chain &chain : model.chains) {
         for (gemmi::Residue &res : chain.residues) {
            for (gemmi::Atom &atom : res.atoms) {
               clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
               clipper::Coord_orth trans_pos = co.transform(rtop);
               atom.pos.x = trans_pos.x();
               atom.pos.y = trans_pos.y();
               atom.pos.z = trans_pos.z();
            }
         }
      }
   }
}

} // namespace coot