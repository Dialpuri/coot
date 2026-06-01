#include "function.hh"
#include <limits>
#include <algorithm>

std::vector<float> set_rainbow_colours_gemmi(const gemmi::Structure& st) {
    std::vector<float> rainbow_vals;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            int nres = static_cast<int>(chain.residues.size());
            int seq_no_max = -std::numeric_limits<int>::max();
            int seq_no_min = std::numeric_limits<int>::max();

            // First pass: find min/max residue indices for standard non-HOH residues
            for (int ires = 0; ires < nres; ires++) {
                const gemmi::Residue& residue = chain.residues[ires];
                std::string res_name = residue.name;
                if (res_name != "HOH") {
                    if (coot::util::is_standard_residue_name(res_name)) {
                        int seq_no_this = ires;
                        if (seq_no_this < seq_no_min) seq_no_min = seq_no_this;
                        if (seq_no_this > seq_no_max) seq_no_max = seq_no_this;
                    }
                }
            }

            if (seq_no_max != -std::numeric_limits<int>::max() &&
                seq_no_min != std::numeric_limits<int>::max()) {

                if (seq_no_min < seq_no_max) {
                    // Second pass: compute rainbow values
                    float range = static_cast<float>(seq_no_max - seq_no_min);
                    for (int ires = 0; ires < nres; ires++) {
                        const gemmi::Residue& residue = chain.residues[ires];
                        std::string res_name = residue.name;
                        if (res_name != "HOH") {
                            if (coot::util::is_standard_residue_name(res_name)) {
                                float chain_pos = static_cast<float>(ires) / range;
                                if (chain_pos < 0) chain_pos = 0;
                                if (chain_pos > 1) chain_pos = 1;

                                for (const gemmi::Atom& atom : residue.atoms) {
                                    (void)atom;
                                    // In MMDB: !atom_p->Het means non-HETATM
                                    // In gemmi: residue.het_flag == 'H' means HETATM record
                                    if (residue.het_flag != 'H') {
                                        rainbow_vals.push_back(chain_pos);
                                    } else {
                                        rainbow_vals.push_back(0.88f);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return rainbow_vals;
}