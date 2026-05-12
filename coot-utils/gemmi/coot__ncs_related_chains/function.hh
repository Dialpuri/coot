#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <map>
#include <string>

namespace coot {

// Port of coot::ncs_related_chains from MMDB to gemmi
// Groups chains that have >70% sequence identity
inline std::vector<std::vector<const gemmi::Chain*>>
ncs_related_chains_gemmi(const gemmi::Model& model) {
    std::vector<std::vector<const gemmi::Chain*>> v;

    auto chains_match = [](const gemmi::Chain* c1, const gemmi::Chain* c2) {
        int n_res_1 = static_cast<int>(c1->residues.size());
        int n_res_2 = static_cast<int>(c2->residues.size());
        
        std::map<int, std::string> m1;
        std::map<int, std::string> m2;
        for (int ires = 0; ires < n_res_1; ires++) {
            const gemmi::Residue& residue = c1->residues[ires];
            std::string res_name = residue.name;
            int res_no = residue.seqid.num.value;
            if (res_name != "HOH")
                m1[res_no] = res_name;
        }
        for (int ires = 0; ires < n_res_2; ires++) {
            const gemmi::Residue& residue = c2->residues[ires];
            std::string res_name = residue.name;
            int res_no = residue.seqid.num.value;
            if (res_name != "HOH")
                m2[res_no] = res_name;
        }

        std::map<int, std::string>::const_iterator it_1;
        unsigned int n_count = 0;
        unsigned int n_match = 0;
        for (it_1 = m1.begin(); it_1 != m1.end(); ++it_1) {
            n_count++;
            std::map<int, std::string>::const_iterator it_2 = m2.find(it_1->first);
            if (it_2 != m2.end()) {
                if (it_2->second == it_1->second)
                    n_match++;
            }
        }
        if (n_count > 0) {
            float f1 = static_cast<float>(n_count);
            float f2 = static_cast<float>(n_match);
            if (f2 / f1 > 0.7)
                return true;
        }
        return false;
    };

    for (const gemmi::Chain& chain : model.chains) {
        if (!chain.residues.empty()) {
            bool found_a_match = false;
            for (size_t iv1 = 0; iv1 < v.size(); ++iv1) {
                std::vector<const gemmi::Chain*>& vv = v[iv1];
                const gemmi::Chain* vv_chain = vv[0];
                bool cm = chains_match(vv_chain, &chain);
                if (cm) {
                    vv.push_back(&chain);
                    found_a_match = true;
                    break;
                }
            }
            if (!found_a_match) {
                v.push_back({&chain});
            }
        }
    }

    return v;
}

} // namespace coot