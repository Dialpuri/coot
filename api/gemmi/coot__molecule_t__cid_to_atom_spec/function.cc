#include "function.hh"
#include <sstream>
#include <algorithm>

namespace coot {

static std::pair<bool, atom_spec_t> parse_cid(const std::string& cid) {
   atom_spec_t spec;
   // CID format: //CHAIN/RESNO/ATOM or //CHAIN/RESNO(ins)/ATOM
   // Example: //A/10/CA, //A/1/CA, //A/2/O, //A/3/CB
   
   if (cid.size() < 5 || cid[0] != '/' || cid[1] != '/') {
      return {false, spec};
   }
   
   // Skip leading //
   size_t pos = 2;
   
   // Parse chain ID
   size_t chain_end = cid.find('/', pos);
   if (chain_end == std::string::npos) return {false, spec};
   spec.chain_id = cid.substr(pos, chain_end - pos);
   if (spec.chain_id.empty()) return {false, spec};
   
   // Parse residue number (and possibly insertion code)
   pos = chain_end + 1;
   size_t res_end = cid.find('/', pos);
   if (res_end == std::string::npos) return {false, spec};
   
   std::string res_str = cid.substr(pos, res_end - pos);
   if (res_str.empty()) return {false, spec};
   
   // Extract residue number and insertion code
   int res_no = 0;
   std::string res_no_str;
   size_t ins_pos = res_str.find_first_of("0123456789");
   if (ins_pos != std::string::npos) {
      res_no_str = res_str.substr(ins_pos);
      spec.ins_code = res_str.substr(0, ins_pos);
      try {
         res_no = std::stoi(res_no_str);
      } catch (...) {
         return {false, spec};
      }
   } else {
      return {false, spec};
   }
   spec.res_no = res_no;
   
   // Parse atom name
   pos = res_end + 1;
   if (pos >= cid.size()) return {false, spec};
   spec.atom_name = cid.substr(pos);
   if (spec.atom_name.empty()) return {false, spec};
   
   // Trim trailing whitespace from atom_name if needed
   while (!spec.atom_name.empty() && spec.atom_name.back() == ' ') {
      spec.atom_name.pop_back();
   }
   
   spec.model_number = 1;
   spec.alt_conf = "";
   
   return {true, spec};
}

static gemmi::Residue* find_residue_by_seqid(const gemmi::Model& model, const std::string& chain_id, int res_no) {
   for (const auto& chain : model.chains) {
      if (chain.name == chain_id) {
         for (const auto& res : chain.residues) {
            if (res.seqid.num.value == res_no) {
               return const_cast<gemmi::Residue*>(&res);
            }
         }
      }
   }
   return nullptr;
}

std::pair<bool, atom_spec_t>
cid_to_atom_spec_gemmi(const gemmi::Structure& st, const std::string& cid) {
   if (st.models.empty()) {
      return {false, atom_spec_t()};
   }
   
   // Parse the CID to extract spec
   auto parsed = parse_cid(cid);
   if (!parsed.first) {
      return {false, atom_spec_t()};
   }
   
   // Find the residue in the structure
   gemmi::Residue* res = find_residue_by_seqid(st.models[0], parsed.second.chain_id, parsed.second.res_no);
   if (!res) {
      return {false, atom_spec_t()};
   }
   
   // Find the atom in the residue
   for (const auto& atom : res->atoms) {
      // Compare atom name - need to handle padding like MMDB does
      std::string atom_name_stripped = atom.name;
      // MMDB stores atom names padded with spaces, gemmi stores them trimmed
      // Format atom name to 4 characters left-padded (like PDB format)
      std::string formatted_name = atom_name_stripped;
      if (formatted_name.size() < 4) {
         formatted_name = std::string(4 - formatted_name.size(), ' ') + formatted_name;
      }
      
      // The CID atom name might be like "CA" (2 chars) or " O  " (4 chars)
      // Let's normalize both for comparison
      std::string cid_atom = parsed.second.atom_name;
      // Remove leading/trailing spaces from CID atom name for comparison
      while (!cid_atom.empty() && cid_atom.front() == ' ') cid_atom.erase(0, 1);
      while (!cid_atom.empty() && cid_atom.back() == ' ') cid_atom.pop_back();
      
      if (atom_name_stripped == cid_atom) {
         atom_spec_t result;
         result.chain_id = parsed.second.chain_id;
         result.res_no = parsed.second.res_no;
         result.ins_code = parsed.second.ins_code;
         // Format atom name to 4 characters for atom_spec_t
         // Pad with spaces on both sides: " CA " not "  CA"
         std::string padded = atom_name_stripped;
         if (padded.size() == 1) {
            padded = " " + padded + "  ";
         } else if (padded.size() == 2) {
            padded = " " + padded + " ";
         } else if (padded.size() == 3) {
            padded = " " + padded;
         }
         result.atom_name = padded;
         result.alt_conf = "";
         result.model_number = 1;
         return {true, result};
      }
   }
   
   return {false, atom_spec_t()};
}

}