#include "function.hh"

namespace coot {
namespace protein_geometry {

bool are_linked_in_order_gemmi(const gemmi::CRA& first_cra,
                               const gemmi::CRA& second_cra,
                               const mmdb::Link* link) {
   bool linked = false;
   std::string link_chain_id_1(link->chainID1);
   std::string link_chain_id_2(link->chainID2);
   std::string chain_id_1 = first_cra.chain->name;
   std::string chain_id_2 = second_cra.chain->name;
   int resno_1 = first_cra.residue->seqid.num.value;
   int resno_2 = second_cra.residue->seqid.num.value;
   if (link_chain_id_1 == chain_id_1) {
      if (link_chain_id_2 == chain_id_2) {
         int link_reso_1 = link->seqNum1;
         int link_reso_2 = link->seqNum2;
         if (link_reso_1 == resno_1) {
            if (link_reso_2 == resno_2) {
               std::string link_ins_code_1(link->insCode1);
               std::string link_ins_code_2(link->insCode2);
               std::string ins_code_1(1, first_cra.residue->seqid.icode);
               std::string ins_code_2(1, second_cra.residue->seqid.icode);
               if (link_ins_code_1 == ins_code_1) {
                  if (link_ins_code_2 == ins_code_2) {
                     linked = true;
                  }
               }
            }
         }
      }
   }
   return linked;
}

} // namespace protein_geometry
} // namespace coot