#include "function.hh"

std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>
coot::util::peptide_C_N_pairs_gemmi(gemmi::Chain* chain_p) {

    double dist_crit = 2.5;

    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> v;
    if (chain_p) {
        int nres = static_cast<int>(chain_p->residues.size());
        for (int ires = 0; ires < (nres - 1); ires++) {
            gemmi::Atom* c_first = nullptr;
            gemmi::Atom* n_next = nullptr;

            gemmi::Residue* residue_p_1 = &chain_p->residues[ires];
            gemmi::Residue* residue_p_2 = &chain_p->residues[ires + 1];

            int seq_num_1 = residue_p_1->seqid.num.value;
            int seq_num_2 = residue_p_2->seqid.num.value;

            // just test for tandem residues, insertion codes will cause a mess
            if (seq_num_2 == (seq_num_1 + 1)) {
                for (gemmi::Atom& at_1 : residue_p_1->atoms) {
                    if (at_1.name == "C") {
                        c_first = &at_1;
                        break;
                    }
                }
                for (gemmi::Atom& at_2 : residue_p_2->atoms) {
                    if (at_2.name == "N") {
                        n_next = &at_2;
                        break;
                    }
                }
                if (c_first) {
                    if (n_next) {
                        // isTer() returns true for "not calculated" atoms.
                        // !isTer() means the atom IS calculated.
                        // In gemmi: calc_flag == NotSet means calculated.
                        if (c_first->calc_flag == gemmi::CalcFlag::NotSet &&
                            n_next->calc_flag == gemmi::CalcFlag::NotSet) {
                            char alt_conf_1 = c_first->altloc;
                            char alt_conf_2 = n_next->altloc;
                            // In MMDB altLoc.empty() means no alt location.
                            // In gemmi altloc == '\0' means no alt location.
                            if (alt_conf_1 == '\0' || alt_conf_2 == '\0' || alt_conf_1 == alt_conf_2) {
                                clipper::Coord_orth pt_1 = coot::co_gemmi(c_first);
                                clipper::Coord_orth pt_2 = coot::co_gemmi(n_next);
                                double d = clipper::Coord_orth::length(pt_1, pt_2);
                                if (d < dist_crit) {
                                    std::string rn_1(residue_p_1->name);
                                    std::string rn_2(residue_p_2->name);
                                    if (is_standard_amino_acid_name(rn_1)) {
                                        if (is_standard_amino_acid_name(rn_2)) {
                                            std::pair<gemmi::Atom*, gemmi::Atom*> p(c_first, n_next);
                                            v.push_back(p);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return v;
}