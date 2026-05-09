#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

inline bool residue_is_nucleic_acid_gemmi(const std::string& cid, const gemmi::Structure& st) {
    bool status = false;
    
    // Parse CID format: //chain/resseq
    // Example: //A/10, //B/20
    std::string chain_id;
    int resseq = 0;
    char icode = ' ';
    
    // Parse CID
    size_t slash1 = cid.find('/');
    if (slash1 != std::string::npos) {
        size_t slash2 = cid.find('/', slash1 + 1);
        if (slash2 != std::string::npos) {
            chain_id = cid.substr(slash1 + 1, slash2 - slash1 - 1);
            size_t slash3 = cid.find('/', slash2 + 1);
            std::string res_part = (slash3 != std::string::npos) ? 
                                   cid.substr(slash2 + 1, slash3 - slash2 - 1) :
                                   cid.substr(slash2 + 1);
            
            // Parse residue sequence number and insertion code
            size_t pos = 0;
            while (pos < res_part.size() && std::isdigit(res_part[pos])) {
                resseq = resseq * 10 + (res_part[pos] - '0');
                pos++;
            }
            if (pos < res_part.size()) {
                icode = res_part[pos];
            }
        }
    }
    
    // Find matching residues in the structure
    std::vector<gemmi::Residue*> matching_residues;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name == chain_id) {
                for (const gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == resseq && res.seqid.icode == icode) {
                        matching_residues.push_back(const_cast<gemmi::Residue*>(&res));
                    }
                }
            }
        }
    }
    
    if (!matching_residues.empty()) {
        bool all_true = true;
        for (gemmi::Residue* res : matching_residues) {
            // Check if residue is a nucleotide
            bool is_nucleotide = false;
            
            // Check using residue name
            std::string res_name = res->name;
            if (res_name == "A" || res_name == "G" || res_name == "C" || 
                res_name == "T" || res_name == "U" || res_name == "DA" || 
                res_name == "DG" || res_name == "DC" || res_name == "DT" || 
                res_name == "DU" || res_name == "RA" || res_name == "RG" || 
                res_name == "RC" || res_name == "RT" || res_name == "RU") {
                is_nucleotide = true;
            }
            
            if (!is_nucleotide) {
                all_true = false;
                break;
            }
        }
        if (all_true) status = true;
    }
    
    return status;
}

}
}