#pragma once
#include <string>
#include <utility>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct residue_spec_t {
   int model_number;
   std::string chain_id;
   int res_no;
   std::string ins_code;
   int int_user_data;
   float float_user_data;
   std::string string_user_data;
   
   residue_spec_t() : model_number(0), res_no(0), int_user_data(0), float_user_data(0) {}
   
   explicit residue_spec_t(int r) : res_no(r), model_number(0), int_user_data(0), float_user_data(0) {}
   
   residue_spec_t(const std::string &chain_in, int r, const std::string &ins_code_in = "")
      : chain_id(chain_in), res_no(r), ins_code(ins_code_in), model_number(0), int_user_data(0), float_user_data(0) {}
};

// Helper function to parse CID string in format //<chain_id>/<res_no>
inline std::pair<std::string, int> parse_cid_residue(const std::string &cid) {
   // Expect format: //<chain_id>/<res_no>
   if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/')
      return std::make_pair("", -1);
   
   // Find the second slash
   size_t pos2 = cid.find('/', 2);
   if (pos2 == std::string::npos)
      return std::make_pair("", -1);
   
   std::string chain_id = cid.substr(2, pos2 - 2);
   if (chain_id.empty())
      return std::make_pair("", -1);
   
   // Parse residue number after the second slash
   std::string num_str = cid.substr(pos2 + 1);
   if (num_str.empty())
      return std::make_pair("", -1);
   
   int res_no = std::stoi(num_str);
   return std::make_pair(chain_id, res_no);
}

// Function to convert CID to residue_spec_t using gemmi
inline std::pair<bool, residue_spec_t> cid_to_residue_spec_gemmi(const gemmi::Structure &st, const std::string &cid) {
   bool status = false;
   residue_spec_t rs;
   
   auto parsed = parse_cid_residue(cid);
   const std::string &chain_id = parsed.first;
   int res_no = parsed.second;
   
   if (chain_id.empty() || res_no < 0)
      return std::make_pair(false, rs);
   
   // Search through models and chains
   for (const auto &model : st.models) {
      for (const auto &chain : model.chains) {
         if (chain.name == chain_id) {
            for (const auto &residue : chain.residues) {
               if (residue.seqid.num.value == res_no) {
                  rs.chain_id = chain_id;
                  rs.res_no = res_no;
                  // Use ins_code only if it's not a space
                  char icode = residue.seqid.icode;
                  rs.ins_code = (icode == ' ') ? "" : std::string(1, icode);
                  status = true;
                  return std::make_pair(status, rs);
               }
            }
         }
      }
   }
   
   return std::make_pair(status, rs);
}

}