#pragma once
#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {

struct atom_overlaps_container_t {
    // Constructor - stub for now
    atom_overlaps_container_t(const gemmi::Model* model, void* geom, bool flag, float dist1, float dist2)
        : model_(model), flag_(flag), dist1_(dist1), dist2_(dist2) {}

    // Ported function - takes CRA instead of bare Residue*
    bool are_bonded_residues_gemmi(const gemmi::CRA& cra1, const gemmi::CRA& cra2) const {
        bool r = false;

        if (cra1.residue && cra2.residue) {
            // Check if they are in the same chain
            if (cra1.chain && cra2.chain && cra1.chain->name == cra2.chain->name) {
                // Check sequence distance
                int seq1 = cra1.residue->seqid.num.value;
                int seq2 = cra2.residue->seqid.num.value;
                
                if (std::abs(seq1 - seq2) < 2) {
                    std::string res_name_1 = cra1.residue->name;
                    std::string res_name_2 = cra2.residue->name;
                    if (res_name_1 != "HOH" && res_name_2 != "HOH") {
                        r = true;
                    }
                }
            }

            if (!r) {
                // Check by serial index (continuous residues in file)
                // We need to find the indices - iterate through all residues
                // to find their positions
                int idx1 = -1, idx2 = -1;
                int current_idx = 0;
                
                for (const auto& chain : model_->chains) {
                    for (const auto& res : chain.residues) {
                        if (&res == cra1.residue) idx1 = current_idx;
                        if (&res == cra2.residue) idx2 = current_idx;
                        current_idx++;
                    }
                }
                
                if (idx1 >= 0 && idx2 >= 0) {
                    if (idx2 > idx1) {
                        if ((idx2 - idx1) == 1) {
                            r = true;
                        }
                    }
                    if (idx1 > idx2) {
                        if ((idx1 - idx2) == 1) {
                            r = true;
                        }
                    }
                }
            }
        }
        return r;
    }

private:
    const gemmi::Model* model_;
    bool flag_;
    float dist1_;
    float dist2_;
};

} // namespace coot