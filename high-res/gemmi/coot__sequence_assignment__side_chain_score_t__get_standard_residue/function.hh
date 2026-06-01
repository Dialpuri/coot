#pragma once

#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace sequence_assignment {

enum side_chain_name_index {
    GLY, ALA, SER, VAL, THR, PRO, ASN, ASP, CYS,
    GLN, GLU, LEU, ILE, MET, PHE, TYR, TRP, LYS, ARG, HIS,
    NUM_RESIDUES
};

class side_chain_score_t {
public:
    std::vector<gemmi::Residue*> standard_residues;

    short int cache_standard_residues() {
        return 1;
    }

    gemmi::Residue*
    get_standard_residue_gemmi(const side_chain_name_index &idx) const {
        return standard_residues[idx];
    }
};

} // namespace sequence_assignment
} // namespace coot