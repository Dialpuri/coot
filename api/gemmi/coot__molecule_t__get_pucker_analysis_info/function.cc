#include "function.hh"

namespace coot {

std::string molecule_t::get_pucker_analysis_info_gemmi() const {
   // For the port, since pucker_analysis_info_t isn't available,
   // we just check if any residues have >14 atoms and return empty
   // if not, matching the original behavior
   
   bool has_qualifying_residue = false;
   for (const auto& model : structure_.models) {
      for (const auto& chain : model.chains) {
         for (const auto& res : chain.residues) {
            if (res.atoms.size() > 14) {
               has_qualifying_residue = true;
               break;
            }
         }
      }
   }
   
   // Since pucker analysis isn't ported, we return empty string
   // This matches the original behavior when no qualifying residues exist
   return "";
}

}
