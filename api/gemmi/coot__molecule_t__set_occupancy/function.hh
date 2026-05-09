#pragma once
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {
namespace molecule_t {

// Set occupancy for atoms matching CID
// CID format: //CHAIN/RESNUM (e.g., "//A/10")
inline void set_occupancy(gemmi::Structure& st, const std::string& cid, float occ_new) {
    if (!st.models.empty()) {
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                // Parse CID to extract chain and residue number
                if (cid.size() >= 4 && cid.substr(0, 2) == "//") {
                    size_t slash_pos = cid.find('/', 2);
                    if (slash_pos != std::string::npos) {
                        std::string chain_id = cid.substr(2, slash_pos - 2);
                        std::string res_str = cid.substr(slash_pos + 1);
                        
                        if (chain.name == chain_id) {
                            // Compare residue sequence number
                            if (std::to_string(residue.seqid.num.value) == res_str) {
                                for (gemmi::Atom& atom : residue.atoms) {
                                    // Skip TERMINATED atoms (named "TER")
                                    if (atom.name != "TER") {
                                        atom.occ = occ_new;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

inline void set_occupancy_gemmi(gemmi::Structure& st, const std::string& cid, float occ_new) {
    set_occupancy(st, cid, occ_new);
}

} // namespace molecule_t
} // namespace coot