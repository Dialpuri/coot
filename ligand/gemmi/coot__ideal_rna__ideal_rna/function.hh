#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace ideal_rna {

gemmi::Structure ideal_rna(const std::string &RNA_or_DNA, const std::string &form,
                           short int single_stranded_flag,
                           const std::string &sequence, gemmi::Structure *standard_residues);

} // namespace ideal_rna
} // namespace coot