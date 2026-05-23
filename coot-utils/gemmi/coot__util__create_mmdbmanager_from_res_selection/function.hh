#pragma once

#include <string>
#include <utility>
#include <vector>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Port of create_mmdbmanager_from_res_selection to gemmi
// Takes a gemmi Structure and vector of Residue pointers, returns new Structure
std::pair<gemmi::Structure, int>
create_mmdbmanager_from_res_selection_gemmi(const gemmi::Structure& st,
                                            const std::vector<gemmi::Residue*>& residues,
                                            int have_flanking_residue_at_start,
                                            int have_flanking_residue_at_end,
                                            const std::string& altconf,
                                            const std::string& chain_id_1,
                                            short int residue_from_alt_conf_split_flag);

}} // namespace coot::util