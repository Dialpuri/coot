#pragma once

#include <gemmi/model.hpp>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__ideal_rna__get_standard_residue_instance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__mutate_base/gemmi/function.hh"

namespace coot {
namespace ideal_rna {

inline int mutate_res_gemmi(
    gemmi::Residue *res,
    char base,
    bool is_dna_flag,
    const gemmi::Model &standard_residues)
{
    // we need to get instances of bases from the standard residues
    int status = 0;

    std::string residue_type = "None";

    if (is_dna_flag) {
        if (base == 'a')
            residue_type = "DA";
        if (base == 'g')
            residue_type = "DG";
        if (base == 't')
            residue_type = "DT";
        if (base == 'c')
            residue_type = "DC";
    } else {
        if (base == 'a')
            residue_type = "A";
        if (base == 'g')
            residue_type = "G";
        if (base == 'u')
            residue_type = "U";
        if (base == 'c')
            residue_type = "C";
    }

    if (residue_type != "None") {
        gemmi::Residue *std_res = get_standard_residue_instance_gemmi(residue_type, standard_residues);
        if (std_res) {
            coot::util::mutate_base_gemmi(*res, *std_res, true, false, 0.0f);
            status = 1;
            delete std_res; // get_standard_residue_instance_gemmi returns a deep-copied pointer
        }
    }
    return status;
}

} // namespace ideal_rna
} // namespace coot