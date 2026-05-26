#pragma once

#include <vector>
#include <set>
#include <queue>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "geometry/residue-and-atom-specs-gemmi.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residues_near_residue__006fec/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__lidia_utils__get_residue/gemmi/function.hh"

namespace coot {

// Convert gemmi Residue to residue_spec_t
inline residue_spec_t residue_to_spec(const gemmi::Residue& res, const gemmi::Chain& chain) {
    std::string chain_id = chain.name;
    int res_no = res.seqid.num.value;
    char icode = res.seqid.icode;
    // Normalize insertion code: gemmi uses ' ' for missing, but residue_spec_t expects " " for missing
    std::string ins_code = (icode == ' ') ? " " : std::string(1, icode);
    return residue_spec_t(chain_id, res_no, ins_code);
}

// Port of coot::simple_residue_tree to gemmi
inline std::vector<gemmi::Residue*>
simple_residue_tree_gemmi(gemmi::Residue* residue_centre, gemmi::Structure& st, float close_dist_max) {
    double dist_crit = close_dist_max;
    std::vector<gemmi::Residue*> v;
    std::set<gemmi::Residue*> s;

    // Get the model and find the chain containing the center residue
    if (st.models.empty()) {
        return v;
    }
    
    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = nullptr;
    for (auto& c : model.chains) {
        for (auto& r : c.residues) {
            if (&r == residue_centre) {
                chain = &c;
                break;
            }
        }
        if (chain) break;
    }
    
    if (!chain) {
        return v;
    }

    std::queue<gemmi::Residue*> q;
    q.push(residue_centre);
    s.insert(residue_centre);

    while (q.size()) {
        gemmi::Residue* test_residue = q.front();
        s.insert(test_residue);
        q.pop();

        // Convert to residue_spec_t for the coot function
        residue_spec_t spec = residue_to_spec(*test_residue, *chain);
        
        // Get nearby residues using the gemmi port
        std::vector<residue_spec_t> specs = residues_near_residue_gemmi(spec, st, dist_crit);
        
        for (const auto& rs : specs) {
            // Find the residue in the structure
            gemmi::Residue* r = lidia_utils::get_residue_gemmi(rs, &st);
            if (r && s.find(r) == s.end()) {
                q.push(r);
                s.insert(r);
            }
        }
    }

    std::set<gemmi::Residue*>::const_iterator its;
    for (its=s.begin(); its!=s.end(); its++)
        v.push_back(*its);

    return v;
}

} // namespace coot