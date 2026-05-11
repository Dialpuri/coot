#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace chi_angles {

struct atom_name_pair {
    std::string atom1;
    std::string atom2;
    atom_name_pair(const std::string& a1, const std::string& a2)
        : atom1(a1), atom2(a2) {}
};

struct atom_index_pair {
    int index1;
    int index2;
    atom_index_pair(int i1, int i2) : index1(i1), index2(i2) {}
};

inline std::vector<atom_index_pair>
get_atom_index_pairs_gemmi(const std::vector<atom_name_pair>& atom_name_pairs,
                           const gemmi::Residue& residue) {
    int i_store_index;
    std::vector<atom_index_pair> index_pairs;

    for (unsigned int ipair = 0; ipair < atom_name_pairs.size(); ipair++) {
        i_store_index = -1;
        for (size_t i = 0; i < residue.atoms.size(); i++) {
            const std::string& atomname = residue.atoms[i].name;
            if (atomname == atom_name_pairs[ipair].atom1) {
                i_store_index = static_cast<int>(i);
            }
        }
        if (i_store_index > -1) {
            for (size_t i2 = 0; i2 < residue.atoms.size(); i2++) {
                const std::string& atomname = residue.atoms[i2].name;
                if (atomname == atom_name_pairs[ipair].atom2) {
                    index_pairs.push_back(atom_index_pair(i_store_index, static_cast<int>(i2)));
                }
            }
        }
    }
    return index_pairs;
}

} // namespace chi_angles
} // namespace coot