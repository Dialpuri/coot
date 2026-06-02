#pragma once

#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"

namespace coot { namespace ideal_rna {

// No gemmi port for this — pure string logic, same as MMDB original.
inline std::string residue_name_from_type(const std::string& orig) {
    // The standard-residues.pdb already contains residues with the original names
    // ("A", "DA", "G", "DG", etc.), so no conversion is needed.
    return orig;
}

inline gemmi::Residue* get_standard_residue_instance_gemmi(
    const std::string& residue_type_in,
    const gemmi::Model& standard_residues)
{
    std::string residue_name = residue_name_from_type(residue_type_in);
    std::cout << "get_standard_residue_instance(): in :" << residue_type_in
              << ": out :" << residue_name << ":" << std::endl;

    // Search through all chains/residues for a matching residue name
    for (const auto& chain : standard_residues.chains) {
        for (const auto& res : chain.residues) {
            if (res.name == residue_name) {
                return coot::util::deep_copy_this_residue_gemmi(&res);
            }
        }
    }

    std::cout << "This should never happen - "
              << "badness in ideal_rna::get_standard_residue_instance(), "
              << "we selected 0 residues looking for residues of type :"
              << residue_name << ": from :" << residue_type_in << ":\n";

    return nullptr;
}

}} // namespace coot::ideal_rna