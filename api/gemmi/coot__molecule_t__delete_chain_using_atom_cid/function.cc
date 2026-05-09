#include "function.hh"
#include <string>
#include <gemmi/pdb.hpp>

namespace coot {

int
delete_chain_using_atom_cid_gemmi(gemmi::Structure& st, const std::string& cid) {
   int done = 0;
   
   // Parse CID to get chain_id
   // CID format: "//CHAIN/RESNO[ICODE]/ATOMNAME"
   std::string chain_id;
   if (cid.size() >= 4 && cid[0] == '/' && cid[1] == '/') {
      size_t pos1 = 2;
      size_t pos2 = cid.find('/', pos1);
      if (pos2 != std::string::npos) {
         chain_id = cid.substr(pos1, pos2 - pos1);
      }
   }
   
   if (!chain_id.empty()) {
      // Only operate on first model (matching MMDB behavior)
      if (!st.models.empty()) {
         gemmi::Model& model = st.models[0];
         
         // Find chain to delete
         for (size_t ichain = 0; ichain < model.chains.size(); ++ichain) {
            if (model.chains[ichain].name == chain_id) {
               model.chains.erase(model.chains.begin() + ichain);
               done = 1;
               break;
            }
         }
      }
   }
   
   return done;
}

}