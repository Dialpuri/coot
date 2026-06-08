#include "function.hh"

namespace coot {
namespace util {

void backrub_residue_triple_t::trim_residue_atoms_generic_gemmi(
    gemmi::Residue* residue_p,
    std::vector<std::string> keep_atom_vector,
    bool use_keep_atom_vector) {

    if (!residue_p) return;

    std::vector<int> delete_atom_index_vec;

    for (int i = 0; i < static_cast<int>(residue_p->atoms.size()); i++) {
        const gemmi::Atom& atom = residue_p->atoms[i];
        std::string atom_name(atom.name);
        // gemmi: altloc is a char (' ' means no alt location, like "" in MMDB)
        std::string atom_alt_conf = (atom.altloc != ' ') ? std::string(1, atom.altloc) : "";

        bool delete_this_atom_flag = true;
        if (use_keep_atom_vector) {
            for (unsigned int ikeep = 0; ikeep < keep_atom_vector.size(); ikeep++) {
                if (atom_name == keep_atom_vector[ikeep]) {
                    delete_this_atom_flag = false;
                    break;
                }
            }
        } else {
            delete_this_atom_flag = false;
        }

        if (delete_this_atom_flag) {
            delete_atom_index_vec.push_back(i);
        } else {
            if ((atom_alt_conf != alt_conf) && (atom_alt_conf != "")) {
                delete_atom_index_vec.push_back(i);
            }
        }
    }

    if (!delete_atom_index_vec.empty()) {
        // Delete in reverse order to avoid index-shifting issues
        for (int i = static_cast<int>(delete_atom_index_vec.size()) - 1; i >= 0; i--) {
            residue_p->atoms.erase(residue_p->atoms.begin() + delete_atom_index_vec[i]);
        }
    }
}

} // namespace util
} // namespace coot
