#include "function.hh"
#include <gemmi/seqid.hpp>

namespace coot {

std::pair<coot::atom_spec_t, coot::atom_spec_t>
link_atoms_gemmi(const gemmi::Connection& link, gemmi::Model* model_p) {

   auto get_res_no = [](const gemmi::ResidueId& rid) -> int {
       return rid.seqid.num.has_value() ? rid.seqid.num.value : 0;
   };
   auto get_ic = [](const gemmi::ResidueId& rid) -> std::string {
       return std::string(1, rid.seqid.icode);
   };

   atom_spec_t a1(link.partner1.chain_name,
                  get_res_no(link.partner1.res_id),
                  get_ic(link.partner1.res_id),
                  link.partner1.atom_name,
                  std::string(1, link.partner1.altloc));
   atom_spec_t a2(link.partner2.chain_name,
                  get_res_no(link.partner2.res_id),
                  get_ic(link.partner2.res_id),
                  link.partner2.atom_name,
                  std::string(1, link.partner2.altloc));

   if (model_p) {
      int mn = model_p->num;
      a1.model_number = mn;
      a2.model_number = mn;
   }

   return {a1, a2};
}

} // namespace coot