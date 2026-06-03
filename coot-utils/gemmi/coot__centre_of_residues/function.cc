#include "function.hh"

namespace coot {

std::pair<bool, gemmi::Vec3>
centre_of_residues_gemmi(const std::vector<gemmi::Residue*> &residues) {

   bool status = false;
   double xs = 0, ys = 0, zs = 0;
   int n_atoms = 0;

   for (auto* res : residues) {
      for (auto& atom : res->atoms) {
         xs += atom.pos.x;
         ys += atom.pos.y;
         zs += atom.pos.z;
         n_atoms++;
      }
   }

   if (n_atoms) {
      status = true;
      double dna = static_cast<double>(n_atoms);
      return {status, gemmi::Vec3(xs / dna, ys / dna, zs / dna)};
   }

   return {status, gemmi::Vec3(0, 0, 0)};

}

} // namespace coot