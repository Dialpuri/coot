#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__deep_copy_this_residue_old_style/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__copy_segid/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_last_residue_in_chain/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__pdbcleanup_serial_residue_numbers/gemmi/function.hh"

// Inlined from coot::find_serial_number_for_insert gemmi port
namespace coot {

inline std::pair<int, gemmi::CRA>
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

} // namespace coot

void insert_coords_gemmi(gemmi::Structure& target, const gemmi::Model& frag_model) {
   bool inserted = false;

   for (const auto& frag_chain : frag_model.chains) {
      for (const auto& frag_res : frag_chain.residues) {
         gemmi::CRA frag_cra{const_cast<gemmi::Chain*>(&frag_chain), const_cast<gemmi::Residue*>(&frag_res), nullptr};

         // Find matching chain in target
         bool chain_matched = false;
         for (auto& target_chain : target.models[0].chains) {
            if (frag_chain.name == target_chain.name) {
               chain_matched = true;

               gemmi::Residue *res = coot::deep_copy_this_residue_old_style_gemmi(frag_cra, "", 1, false);

               if (res) {
                  auto serial_number = coot::find_serial_number_for_insert_gemmi(
                     frag_res.seqid.num.value,
                     std::string(1, frag_res.seqid.icode),
                     target_chain.name,
                     target.models[0]);

                  if (serial_number.first != -1) {
                     auto it = target_chain.residues.begin() + serial_number.first;
                     target_chain.residues.insert(it, *res);
                     if (serial_number.second.residue) {
                        coot::copy_segid_gemmi(*serial_number.second.residue, target_chain.residues[serial_number.first]);
                     }
                     inserted = true;
                  } else {
                     auto last_residue = coot::util::get_last_residue_in_chain_gemmi(target_chain);
                     target_chain.residues.push_back(*res);
                     if (last_residue) {
                        coot::copy_segid_gemmi(*last_residue, target_chain.residues.back());
                     }
                     inserted = true;
                  }
                  delete res;
               }
               break;
            }
         }

         if (!chain_matched) {
            // Create a new chain in the target molecule
            gemmi::Chain new_chain;
            new_chain.name = frag_chain.name;

            gemmi::Residue *res = coot::deep_copy_this_residue_old_style_gemmi(frag_cra, "", 1, false);
            if (res) {
               new_chain.residues.push_back(*res);
               target.models[0].chains.push_back(new_chain);
               delete res;
               inserted = true;
            }
         }
      }
   }

   coot::util::pdbcleanup_serial_residue_numbers_gemmi(target);
}