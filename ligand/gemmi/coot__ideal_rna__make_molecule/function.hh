#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <memory>
#include <vector>
#include "clipper/core/coords.h"

namespace coot {
namespace ideal_rna {

// Helper functions (no _gemmi port yet)
void delete_o2_prime_gemmi(gemmi::Residue* res);
bool is_valid_base_gemmi(char c);
clipper::RTop_orth n_turns_gemmi(unsigned int iseq, int len, int form_flag);
int mutate_res_gemmi(gemmi::Residue* res, char base, bool is_dna_flag);
char antisense_base_gemmi(char c, bool is_dna_flag);

// Main function
std::unique_ptr<gemmi::Structure> make_molecule_gemmi(
    const std::string& seq,
    const std::string& RNA_or_DNA,
    const std::string& form,
    bool single_stranded,
    const gemmi::Model& standard_residues);

} // namespace ideal_rna
} // namespace coot