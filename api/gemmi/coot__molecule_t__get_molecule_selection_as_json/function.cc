#include "function.hh"
#include <sstream>
#include <algorithm>

namespace coot {

bool parse_cid(const std::string &cid, std::string &chain_id, int &res_seqnum) {
    // Format: "//<ChainID>/<ResidueNumber>"
    if (cid.size() < 5 || cid[0] != '/' || cid[1] != '/' || cid[3] != '/')
        return false;
    
    chain_id = cid.substr(2, cid.size() - cid.find('/', 2) - 2);
    std::string num_str = cid.substr(cid.find('/', 2) + 1);
    
    try {
        res_seqnum = std::stoi(num_str);
    } catch (...) {
        return false;
    }
    return true;
}

static std::string remove_whitespace(const std::string &s) {
    std::string result;
    result.reserve(s.size());
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c)))
            result += c;
    }
    return result;
}

static std::string escape_json_string(const std::string &s) {
    std::string result = "\"";
    for (char c : s) {
        if (c == '"') result += "\\\"";
        else if (c == '\\') result += "\\\\";
        else result += c;
    }
    result += "\"";
    return result;
}

std::string get_molecule_selection_as_json_gemmi(const gemmi::Structure &st, const std::string &cid) {
    std::string chain_id;
    int res_seqnum;
    if (!parse_cid(cid, chain_id, res_seqnum))
        return "";
    
    // Find the matching residue
    bool found = false;
    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name != chain_id)
                continue;
            for (auto &res : chain.residues) {
                if (res.seqid.num.value == res_seqnum) {
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
        if (found) break;
    }
    
    if (!found)
        return "";
    
    // Build JSON manually
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"models\": [\n";
    
    bool first_model = true;
    for (auto &model : st.models) {
        if (!first_model)
            oss << "    ,\n";
        first_model = false;
        
        oss << "    {\n";
        oss << "      \"chains\": [\n";
        
        bool first_chain = true;
        for (auto &chain : model.chains) {
            if (!first_chain)
                oss << "        ,\n";
            first_chain = false;
            
            oss << "        {\n";
            oss << "          \"chain_id\": " << escape_json_string(chain.name) << ",\n";
            oss << "          \"residues\": [\n";
            
            bool first_res = true;
            for (auto &res : chain.residues) {
                // Only include residues that have atoms matching our selection
                int atom_count = 0;
                for (auto &atom : res.atoms) {
                    if (!atom.element.is_hydrogen()) {
                        atom_count++;
                    }
                }
                
                if (atom_count == 0)
                    continue;
                
                if (!first_res)
                    oss << "          ,\n";
                first_res = false;
                
                oss << "          {\n";
                oss << "            \"name\": " << escape_json_string(res.name) << ",\n";
                oss << "            \"ins_code\": " << escape_json_string(std::string(1, res.seqid.icode)) << ",\n";
                oss << "            \"seqnum\": " << res.seqid.num.value << ",\n";
                oss << "            \"atoms\": [\n";
                
                bool first_atom = true;
                for (auto &atom : res.atoms) {
                    if (!first_atom)
                        oss << "              ,\n";
                    first_atom = false;
                    
                    std::string elem = atom.element.name();
                    std::string name = remove_whitespace(atom.name);
                    
                    oss << "              {\n";
                    oss << "                \"x\": " << atom.pos.x << ",\n";
                    oss << "                \"y\": " << atom.pos.y << ",\n";
                    oss << "                \"z\": " << atom.pos.z << ",\n";
                    oss << "                \"tempFactor\": " << atom.b_iso << ",\n";
                    oss << "                \"occupancy\": " << atom.occ << ",\n";
                    oss << "                \"name\": " << escape_json_string(name) << ",\n";
                    oss << "                \"element\": " << escape_json_string(elem) << "\n";
                    oss << "              }";
                }
                
                oss << "\n            ]\n";
                oss << "          }";
            }
            
            oss << "\n          ]\n";
            oss << "        }";
        }
        
        oss << "\n      ]\n";
        oss << "    }";
    }
    
    oss << "\n  ]\n";
    oss << "}\n";
    
    return oss.str();
}

} // namespace coot
