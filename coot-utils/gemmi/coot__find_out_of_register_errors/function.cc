#include "function.hh"
#include <gemmi/polyheur.hpp>
#include <iostream>

namespace coot {

void find_out_of_register_errors_gemmi(gemmi::Structure* mutations_st, gemmi::Structure* ref_st) {
    // Add hydrogen atoms if missing and set up entities
    gemmi::setup_entities(*mutations_st);
    gemmi::setup_entities(*ref_st);
    
    // We have done mutations post alignment. There might be out-of-register errors
    // (probably are). I want to find them.
    // The ref_st is the reference/deposited structure.
    // The mutated molecule will not have the correct residue numbering
    // and will have mis-built loops (too many or too few residues).
    //
    // So find 5-residue fragments in the reference, and try to find that
    // sequence in the mutations mol. Then find the difference of the positions
    // based around the central residue (say).
    //
    // Which chains? The first chains in both structures

    if (!mutations_st || mutations_st->models.empty()) return;
    if (!ref_st || ref_st->models.empty()) return;

    class sequence_info_t {
    public:
        sequence_info_t(gemmi::Residue* r, int rn, const std::string& resname)
            : residue_p(r), res_no(rn), res_name(resname) {}
        gemmi::Residue* residue_p;
        int res_no;
        std::string res_name;
    };

    auto fill_sequence = [](gemmi::Structure* st) {
        std::vector<sequence_info_t> sequence;
        if (!st->models.empty()) {
            gemmi::Model& model = st->models[0];
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& residue : chain.residues) {
                    int res_no = residue.seqid.num.value;
                    std::string rn = residue.name;
                    sequence_info_t si(&residue, res_no, rn);
                    sequence.push_back(si);
                }
                if (!sequence.empty())
                    break;
            }
        }
        return sequence;
    };

    // return empty vector on failure
    auto make_a_5_residue_fragment = [](int i, const std::vector<sequence_info_t>& si) {
        std::vector<sequence_info_t> s;
        int si_size = si.size();
        int i_start = i - 2;
        if (i_start >= 0) {
            int i_last = i + 2;   // last index in the fragment
            int i_end = i_last + 1;  // first invalid index
            if (i_end < si_size) {
                std::vector<sequence_info_t>::const_iterator it_begin = si.begin() + i_start;
                std::vector<sequence_info_t>::const_iterator it_end = si.begin() + i_end;
                s = std::vector<sequence_info_t>(it_begin, it_end);
            }
        }
        return s;
    };
    
    // return empty vector on failure
    auto find_fragment_in_sequence = [](const std::vector<sequence_info_t>& search_frag,
                                        const std::vector<sequence_info_t>& sequence) {
        std::vector<sequence_info_t> s;
        for (size_t i = 2; i < sequence.size(); i++) {
            size_t n_matches = 0;
            if (i < (sequence.size() - 2)) {
                for (size_t j = 0; j < search_frag.size(); j++) {
                    if (search_frag[j].res_name == sequence[i + j - 2].res_name)
                        n_matches++;
                }
            }
            if (n_matches == 5) {
                if (s.empty()) {
                    std::vector<sequence_info_t>::const_iterator it_begin = sequence.begin() + i - 2;
                    std::vector<sequence_info_t>::const_iterator it_end = sequence.begin() + i + 3;
                    s = std::vector<sequence_info_t>(it_begin, it_end);
                } else {
                    std::cout << "Double hit for fragment sequence returning first hit " << std::endl;
                }
            }
        }
        return s;
    };

    auto get_positions_difference = [](const std::vector<sequence_info_t>& frag_1,
                                       const std::vector<sequence_info_t>& frag_2) {
        float d = 0.0f;
        size_t n_hits = 0;
        for (size_t i = 0; i < frag_1.size(); i++) {
            gemmi::Residue* res = frag_1[i].residue_p;
            gemmi::Atom* at_1 = nullptr;
            for (gemmi::Atom& atom : res->atoms) {
                if (atom.name == "CA") {
                    at_1 = &atom;
                    break;
                }
            }
            if (at_1) {
                res = frag_2[i].residue_p;
                gemmi::Atom* at_2 = nullptr;
                for (gemmi::Atom& atom : res->atoms) {
                    if (atom.name == "CA") {
                        at_2 = &atom;
                        break;
                    }
                }
                if (at_2) {
                    clipper::Coord_orth pos_1 = co_gemmi(at_1);
                    clipper::Coord_orth pos_2 = co_gemmi(at_2);
                    double dd = (pos_2 - pos_1).lengthsq();
                    d += sqrtf(dd);
                    std::cout << "adding d " << sqrtf(dd) << std::endl;
                    n_hits++;
                }
            }
        }
        if (n_hits != 5)
            return static_cast<float>(-1.0);
        return d;
    };

    std::vector<sequence_info_t> ref_sequence = fill_sequence(ref_st);
    std::vector<sequence_info_t> mut_sequence = fill_sequence(mutations_st);

    if (ref_sequence.empty()) {
        std::cout << "empty ref sequence " << std::endl;
        return;
    }
    if (mut_sequence.empty()) {
        std::cout << "empty mut sequence " << std::endl;
        return;
    }

    std::vector<std::pair<gemmi::Residue*, float>> position_difference_results;
    for (size_t i = 2; i < ref_sequence.size(); i++) {
        std::vector<sequence_info_t> five_res_frag = make_a_5_residue_fragment(i, ref_sequence);
        if (five_res_frag.size() == 5) {
            std::vector<sequence_info_t> related_frag = find_fragment_in_sequence(five_res_frag, mut_sequence);
            if (related_frag.size() == 5) {
                gemmi::Residue* r = five_res_frag[2].residue_p;
                std::cout << "Checking differences for residue " << r->seqid.num.value << std::endl;
                float positions_difference = get_positions_difference(five_res_frag, related_frag);
                float d = positions_difference;
                position_difference_results.push_back(std::make_pair(r, d));
            }
        }
    }

    for (size_t i = 0; i < position_difference_results.size(); i++) {
        std::cout << i << "   " << position_difference_results[i].first->seqid.num.value
                  << " " << position_difference_results[i].second << std::endl;
    }
}

} // namespace coot