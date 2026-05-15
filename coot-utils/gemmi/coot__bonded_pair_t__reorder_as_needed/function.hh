#pragma once

#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>

namespace coot {

// Check if a residue is an amino acid (polypeptide).
// gemmi group_idx: 0=standard AA, 1=modified AA.
inline bool is_aminoacid(const gemmi::Residue* r) {
    return r && r->group_idx <= 1;
}

// Check if a residue is a nucleotide (DNA/RNA).
// gemmi group_idx: 2=DNA, 3=RNA, 4=DNA/RNA hybrid.
inline bool is_nucleotide(const gemmi::Residue* r) {
    return r && r->group_idx >= 2 && r->group_idx <= 4;
}

struct bonded_pair_t {
    gemmi::CRA res_1 = {nullptr, nullptr, nullptr};
    gemmi::CRA res_2 = {nullptr, nullptr, nullptr};
    bool is_fixed_first = false;
    bool is_fixed_second = false;
    std::string bond_type;

    bonded_pair_t(gemmi::CRA r1, gemmi::CRA r2, bool f1, bool f2,
                  const std::string& bt)
        : res_1(r1), res_2(r2), is_fixed_first(f1), is_fixed_second(f2),
          bond_type(bt) {}

    void reorder_as_needed_gemmi() {
        if (res_2.residue && res_1.residue &&
            res_2.residue->seqid.num.value < res_1.residue->seqid.num.value) {
            std::string chain_id_1, chain_id_2;
            if (res_1.chain) chain_id_1 = res_1.chain->name;
            if (res_2.chain) chain_id_2 = res_2.chain->name;
            if (chain_id_1 == chain_id_2) {
                if (is_aminoacid(res_1.residue)) {
                    if (is_aminoacid(res_2.residue)) {
                        std::swap(res_1, res_2);
                        std::swap(is_fixed_first, is_fixed_second);
                    }
                }
                if (is_nucleotide(res_1.residue)) {
                    if (is_nucleotide(res_2.residue)) {
                        std::swap(res_1, res_2);
                        std::swap(is_fixed_first, is_fixed_second);
                    }
                }
            }
        }
    }
};

} // namespace coot