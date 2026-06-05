#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <iostream>

namespace coot {
namespace rotamer {

// Trim whitespace from both ends of a string (for PDB atom name normalization)
static inline std::string trim_spaces(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

std::vector<std::vector<int>> rotamer_atom_names_to_indices_gemmi(
    const std::vector<std::vector<std::string>>& residue_rotamer_atoms,
    gemmi::CRA cra)
{
    bool verbose = false;
    std::vector<std::vector<int>> r;

    if (!cra.residue)
        return r;

    int n_residue_atoms = static_cast<int>(cra.residue->atoms.size());

    if (n_residue_atoms > 0) {
        std::vector<std::string> atom_indices(n_residue_atoms);
        for (int iat = 0; iat < n_residue_atoms; iat++) {
            atom_indices[iat] = trim_spaces(cra.residue->atoms[iat].name);
        }

        int ithis_atom;
        // for each chi atom name set:
        for (unsigned int ich = 0; ich < residue_rotamer_atoms.size(); ich++) {
            std::vector<int> single;
            // for each atom name (in the chi set):
            for (unsigned int iat = 0; iat < residue_rotamer_atoms[ich].size(); iat++) {
                ithis_atom = -1;
                // Can we find in atom_indices the atom with the same name as
                // the (dictionary) chi atom name?
                for (int irat = 0; irat < n_residue_atoms; irat++) {
                    if (atom_indices[irat] == residue_rotamer_atoms[ich][iat]) {
                        ithis_atom = irat;
                        break;
                    }
                }
                if (ithis_atom != -1) {
                    single.push_back(ithis_atom);
                }
            }
            if (single.size() == 4)
                r.push_back(single);
            else if (verbose) {
                std::string res_name = cra.residue ? cra.residue->name : "???";
                int seqnum = 0;
                if (cra.residue)
                    seqnum = cra.residue->seqid.num.value;
                std::string chain_id = (cra.chain) ? cra.chain->name : "?";
                std::cout << "PROBLEM in coordinates file? failed to find all atoms in "
                          << "ich number " << ich << " "
                          << res_name << " "
                          << seqnum << " "
                          << chain_id << std::endl;
            }
        }
    }
    return r;
}

} // namespace rotamer
} // namespace coot