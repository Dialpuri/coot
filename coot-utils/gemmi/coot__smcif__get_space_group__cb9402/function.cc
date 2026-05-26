#include "function.hh"
#include <clipper/core/spacegroup.h>

namespace coot {
namespace smcif {

std::pair<bool, clipper::Spacegroup>
get_space_group_gemmi(const std::vector<std::string> &symm_strings) {
   bool status = false;
   std::string symmetry_ops;
   for (unsigned int isym=0; isym<symm_strings.size(); isym++) { 
      symmetry_ops += symm_strings[isym];
      symmetry_ops += " ; ";
   }
   clipper::Spacegroup space_group;
   clipper::Spgr_descr spg_descr(symmetry_ops, clipper::Spgr_descr::Symops);

   if (spg_descr.spacegroup_number() == 0) {
      // Failed.
      
   } else {
      // Happy path
      space_group.init(spg_descr);
      status = true;
   }
   return std::pair<bool,clipper::Spacegroup>(status, space_group);
}

} // namespace smcif
} // namespace coot