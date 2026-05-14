#include "function.hh"
#include <gemmi/chemcomp.hpp>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__deep_copy_this_residue_old_style/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__find_serial_number_for_insert/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__copy_segid/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_last_residue_in_chain/gemmi/function.hh"

namespace coot {

bool insert_coords_internal_gemmi(const gemmi::Model& src_model, gemmi::Model& target_model) {
   bool inserted = false;
   
   for (auto& src_chain : const_cast<gemmi::Model&>(src_model).chains) {
      bool chain_found = false;
      for (auto& target_chain : target_model.chains) {
         if (src_chain.name == target_chain.name) {
            chain_found = true;
            
            for (auto& src_residue : src_chain.residues) {
               gemmi::CRA src_cra{&src_chain, &src_residue, nullptr};
               
               std::pair<int, gemmi::CRA> serial_result = 
                  find_serial_number_for_insert_gemmi(
                     src_residue.seqid.num.value,
                     std::string(1, src_residue.seqid.icode),
                     src_chain.name,
                     target_model
                  );
               
               if (serial_result.first != -1) {
                  gemmi::Residue* new_res_ptr = 
                     deep_copy_this_residue_old_style(src_cra, "", 1, false);
                  
                  if (new_res_ptr) {
                     gemmi::Residue new_res = *new_res_ptr;
                     delete new_res_ptr;
                     
                     auto it = target_chain.residues.begin() + serial_result.first;
                     target_chain.residues.insert(it, new_res);
                     
                     if (serial_result.second.residue) {
                        copy_segid_gemmi(*serial_result.second.residue, new_res);
                     }
                     inserted = true;
                  }
               } else {
                  gemmi::Residue* new_res_ptr = 
                     deep_copy_this_residue_old_style(src_cra, "", 1, false);
                  
                  if (new_res_ptr) {
                     gemmi::Residue new_res = *new_res_ptr;
                     delete new_res_ptr;
                     
                     const gemmi::Residue* last_res = util::get_last_residue_in_chain_gemmi(target_chain);
                     if (last_res) {
                        target_chain.residues.push_back(new_res);
                        copy_segid_gemmi(*last_res, new_res);
                     } else {
                        target_chain.residues.push_back(new_res);
                     }
                     inserted = true;
                  }
               }
            }
            break;
         }
      }
      
      if (!chain_found) {
         gemmi::Chain new_chain;
         new_chain.name = src_chain.name;
         
         for (auto& src_res : src_chain.residues) {
            gemmi::CRA src_cra{&src_chain, &src_res, nullptr};
            gemmi::Residue* new_res_ptr = 
               deep_copy_this_residue_old_style(src_cra, "", 1, false);
            
            if (new_res_ptr) {
               gemmi::Residue new_res = *new_res_ptr;
               delete new_res_ptr;
               new_chain.residues.push_back(new_res);
            }
         }
         
         target_model.chains.push_back(new_chain);
         inserted = true;
      }
   }
   
   return inserted;
}

} // namespace coot