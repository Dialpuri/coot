#include "function.hh"
#include <gemmi/resinfo.hpp>

namespace coot { namespace util {

bool is_standard_amino_acid_name(const std::string& residue_name) {
    gemmi::ResidueInfo& ri = gemmi::find_tabulated_residue(residue_name);
    return ri.is_amino_acid() && ri.is_standard();
}

std::string model_sequence_gemmi(const std::vector<std::pair<gemmi::Residue*, int>>& sa,
                                  bool allow_ligands) {
    std::string s;
    for (size_t i = 0; i < sa.size(); i++) {
        std::string this_residue = "X";
        const std::string& res_name = sa[i].first->name;
        if (is_standard_amino_acid_name(res_name) || allow_ligands) {
            gemmi::ResidueInfo& ri = gemmi::find_tabulated_residue(res_name);
            if (ri.found() && ri.one_letter_code != ' ') {
                this_residue = std::string(1, ri.one_letter_code);
            } else {
                this_residue = "X";
            }
            s += this_residue;
        }
    }
    return s;
}

}} // namespace coot::util