#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

class fragment_container_t {
public:
    class fragment_range_t {
    public:
        std::vector<gemmi::CRA> residues; // ordered residues in range
        std::string chain_id;
        struct residue_spec_t {
            std::string chain_id;
            int res_no;
            char ins_code;
            residue_spec_t(const std::string &chain_id_in, int res_no_in, char ins_code_in)
                : chain_id(chain_id_in), res_no(res_no_in), ins_code(ins_code_in) {}
        };
        residue_spec_t start_res;
        residue_spec_t end_res;
        fragment_range_t(const std::string &chain_id_in,
                         const residue_spec_t &r1, const residue_spec_t &r2)
            : chain_id(chain_id_in), start_res(r1), end_res(r2) {}
    };

    std::vector<fragment_range_t> ranges;
    fragment_container_t() {}
    void add(const fragment_range_t &r) {
        ranges.push_back(r);
    }
};

// "fragment" here means whole chains, or parts of chains if there are chain breaks.
inline fragment_container_t make_fragments_gemmi(const gemmi::Structure &st) {
    fragment_container_t fc;

    if (st.models.empty())
        return fc;

    const gemmi::Model &model = st.models[0];
    
    for (const gemmi::Chain &chain : model.chains) {
        std::string chain_id = chain.name;
        const auto &residues = chain.residues;
        
        if (residues.empty())
            continue;
        
        gemmi::Residue *residue_start = nullptr;
        gemmi::Residue *residue_prev = nullptr;
        std::vector<gemmi::CRA> residues_running;
        
        for (size_t i = 0; i < residues.size(); ++i) {
            gemmi::Residue &res = const_cast<gemmi::Residue&>(residues[i]);
            if (!residue_start) {
                residue_start = &res;
                residues_running.push_back({const_cast<gemmi::Chain*>(&chain), &res, nullptr});
            } else {
                if (residue_prev) {
                    int resno_prev = residue_prev->seqid.num.value;
                    int resno_this = res.seqid.num.value;
                    if (resno_prev != (resno_this - 1)) {
                        fragment_container_t::fragment_range_t fr(
                            chain_id,
                            fragment_container_t::fragment_range_t::residue_spec_t(
                                chain_id, residue_start->seqid.num.value, residue_start->seqid.icode),
                            fragment_container_t::fragment_range_t::residue_spec_t(
                                chain_id, residue_prev->seqid.num.value, residue_prev->seqid.icode)
                        );
                        fr.residues = residues_running;
                        residue_start = &res;
                        fc.add(fr);
                        residues_running.clear();
                    } else {
                        residues_running.push_back({const_cast<gemmi::Chain*>(&chain), &res, nullptr});
                    }
                }
            }
            residue_prev = &res;
        }

        if (!residues_running.empty()) {
            if (residue_start && residue_prev) {
                if (residue_start != residue_prev) {
                    fragment_container_t::fragment_range_t fr(
                        chain_id,
                        fragment_container_t::fragment_range_t::residue_spec_t(
                            chain_id, residue_start->seqid.num.value, residue_start->seqid.icode),
                        fragment_container_t::fragment_range_t::residue_spec_t(
                            chain_id, residue_prev->seqid.num.value, residue_prev->seqid.icode)
                    );
                    fr.residues = residues_running;
                    fc.add(fr);
                }
            }
        }
    }

    return fc;
}

} // namespace coot