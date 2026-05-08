#include "function.hh"

namespace coot {

std::pair<int, gemmi::CRA>
find_serial_number_for_insert_gemmi(int seqnum_for_new,
                                    const std::string& ins_code_for_new,
                                    const std::string& chain_id,
                                    const gemmi::Model& model) {
   int iserial_no = -1;
   gemmi::CRA result_cra{nullptr, nullptr, nullptr};

   for (const gemmi::Chain& chain : model.chains) {
      if (chain.name != chain_id) {
         continue;
      }

      for (size_t ires = 0; ires < chain.residues.size(); ++ires) {
         const gemmi::Residue& residue = chain.residues[ires];
         int diff = residue.seqid.num.value - seqnum_for_new;

         if (diff > 0) {
            iserial_no = static_cast<int>(ires);
            result_cra = {const_cast<gemmi::Chain*>(&chain), 
                          const_cast<gemmi::Residue*>(&residue), 
                          nullptr};
            break;
         } else if (diff == 0) {
            // Only compare insertion codes if ins_code_for_new is not empty
            if (!ins_code_for_new.empty()) {
               char ins_code_this = residue.seqid.icode;
               if (ins_code_this > ins_code_for_new[0]) {
                  iserial_no = static_cast<int>(ires);
                  result_cra = {const_cast<gemmi::Chain*>(&chain), 
                                const_cast<gemmi::Residue*>(&residue), 
                                nullptr};
                  break;
               }
            }
         }
      }
   }
   return std::make_pair(iserial_no, result_cra);
}

}