#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>

namespace coot {
namespace sequence_assignment {

struct high_res_residue_range_t {
    std::string chain_id;
    int res_no_start;
    int res_no_end;
    char ins_code_start;
    char ins_code_end;
};

inline void generate_scores_gemmi(
    gemmi::Model& mol_in,
    const clipper::Xmap<float>& xmap,
    std::map<std::string, std::vector<std::pair<int, char>>>& seq
) {
    // Build a set of assigned residues: chain_id -> set of (res_no, ins_code)
    std::map<std::string, std::set<std::pair<int, char>>> assigned_set;
    for (auto& chain : mol_in.chains) {
        std::string chain_id = chain.name;
        auto seq_it = seq.find(chain_id);
        if (seq_it != seq.end()) {
            for (auto& seq_res : seq_it->second) {
                assigned_set[chain_id].insert(seq_res);
            }
        }
    }

    // ----- find_unassigned_regions logic -----
    float pr_cut = 0.1;
    std::vector<high_res_residue_range_t> urv;
    for (auto& chain : mol_in.chains) {
        std::string chain_id = chain.name;
        int start = -1;
        int end = -1;
        char start_ic = ' ';
        char end_ic = ' ';
        int unassigned_count = 0;
        int total_count = 0;

        for (auto& res : chain.residues) {
            total_count++;
            int res_no = res.seqid.num.value;
            char ins_code = res.seqid.icode;
            bool is_assigned = false;
            auto chain_it = assigned_set.find(chain_id);
            if (chain_it != assigned_set.end()) {
                if (chain_it->second.count(std::make_pair(res_no, ins_code))) {
                    is_assigned = true;
                }
            }

            if (!is_assigned) {
                if (start == -1) {
                    start = res_no;
                    start_ic = ins_code;
                }
                end = res_no;
                end_ic = ins_code;
                unassigned_count++;
            } else {
                if (start != -1 && total_count >= 3 && ((float)unassigned_count / total_count) > pr_cut) {
                    high_res_residue_range_t rr;
                    rr.chain_id = chain_id;
                    rr.res_no_start = start;
                    rr.res_no_end = end;
                    rr.ins_code_start = start_ic;
                    rr.ins_code_end = end_ic;
                    urv.push_back(rr);
                }
                start = -1;
                end = -1;
                unassigned_count = 0;
                total_count = 0;
            }
        }
        // Handle trailing unassigned region
        if (start != -1 && total_count >= 3 && ((float)unassigned_count / total_count) > pr_cut) {
            high_res_residue_range_t rr;
            rr.chain_id = chain_id;
            rr.res_no_start = start;
            rr.res_no_end = end;
            rr.ins_code_start = start_ic;
            rr.ins_code_end = end_ic;
            urv.push_back(rr);
        }
    }

    std::cout << "There were " << urv.size() << " unassigned regions\n";
}

} // namespace sequence_assignment
} // namespace coot