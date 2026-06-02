#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"

namespace coot {
namespace sequence_assignment {

enum side_chain_name_index {
    GLY, ALA, SER, VAL, THR, PRO, ASN, ASP, CYS,
    GLN, GLU, LEU, ILE, MET, PHE, TYR, TRP, LYS, ARG, HIS,
    NUM_RESIDUES
};

inline std::string side_chain_name_index_to_name(side_chain_name_index idx) {
    static const char* names[] = {
        "GLY", "ALA", "SER", "VAL", "THR", "PRO", "ASN", "ASP", "CYS",
        "GLN", "GLU", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "LYS", "ARG", "HIS"
    };
    int i = static_cast<int>(idx);
    if (i >= 0 && i < 20) return names[i];
    return "UNK";
}

class side_chain_score_t {
public:
    std::vector<gemmi::Residue*> standard_residues;

    short int cache_standard_residues_gemmi(gemmi::Structure* structure = nullptr);
};

inline short int side_chain_score_t::cache_standard_residues_gemmi(gemmi::Structure* structure) {
    if (structure == nullptr) {
        return 1;
    }

    standard_residues.resize(20);
    for (int i = 0; i < 20; i++) {
        std::string residue_type = side_chain_name_index_to_name(static_cast<side_chain_name_index>(i));
        gemmi::Residue* found = nullptr;
        for (gemmi::Model& model : structure->models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.name == residue_type) {
                        found = &res;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (found) {
            standard_residues[i] = coot::util::deep_copy_this_residue_gemmi(found);
        }
    }
    return 1;
}

} // namespace sequence_assignment
} // namespace coot