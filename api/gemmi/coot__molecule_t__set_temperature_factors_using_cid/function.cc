#include "function.hh"
#include <sstream>
#include <cstdlib>

namespace coot {
namespace molecule_t {

// Helper function to parse CID string like "//A/1"
// Returns {chain_name, residue_number} or empty on failure
static std::pair<std::string, int> parse_cid(const std::string& cid) {
    // Format: //chain/resnum (e.g., "//A/1", "//B/10")
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/')
        return {"", -1};
    
    // Find the slash after chain ID
    size_t slash_pos = cid.find('/', 2);
    if (slash_pos == std::string::npos || slash_pos == 2)
        return {"", -1};
    
    std::string chain_id = cid.substr(2, slash_pos - 2);
    
    // Parse residue number
    std::string res_num_str = cid.substr(slash_pos + 1);
    if (res_num_str.empty())
        return {"", -1};
    
    int res_num = std::atoi(res_num_str.c_str());
    
    return {chain_id, res_num};
}

void set_temperature_factors_using_cid_gemmi(gemmi::Structure& st,
                                             const std::string& cid,
                                             float temp_fact) {
    auto parsed = parse_cid(cid);
    std::string chain_name = parsed.first;
    int res_num = parsed.second;
    
    if (chain_name.empty() || res_num < 0)
        return;  // Invalid CID, silently do nothing
    
    // Traverse all models, chains, residues to find matching atoms
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_name)
                continue;
            
            for (auto& res : chain.residues) {
                if (res.seqid.num.value != res_num)
                    continue;
                
                // Found matching residue, update all atom tempFactors
                for (auto& atom : res.atoms) {
                    atom.b_iso = temp_fact;
                }
            }
        }
    }
}

}
}