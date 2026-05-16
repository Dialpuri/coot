#include "function.hh"

namespace coot {

std::vector<gemmi::Residue *>
cid_to_residues_gemmi(const std::string &atom_selection_cids, gemmi::Structure &st) {
    std::vector<gemmi::Residue *> v;

    std::set<gemmi::Residue *> residue_set;
    std::vector<std::string> cid_v = coot::split_string(atom_selection_cids, "||");

    if (!cid_v.empty()) {
        for (const auto &cid : cid_v) {
            // Parse CID format: //chain/residuenumber
            if (cid.size() > 2 && cid[0] == '/' && cid[1] == '/') {
                std::string chain_id = cid.substr(2, cid.find('/', 2) - 2);
                std::string resid_str = cid.substr(cid.find('/', 2) + 1);
                try {
                    int resid = std::stoi(resid_str);
                    for (auto &model : st.models) {
                        for (auto &chain : model.chains) {
                            if (chain.name == chain_id) {
                                for (auto &res : chain.residues) {
                                    if (res.seqid.num.value == resid) {
                                        residue_set.insert(&res);
                                    }
                                }
                            }
                        }
                    }
                } catch (...) {
                    // Invalid number, skip
                }
            }
        }
    }

    v.reserve(residue_set.size());
    for (auto it = residue_set.begin(); it != residue_set.end(); ++it) {
        v.push_back(*it);
    }

    return v;
}

}  // namespace coot