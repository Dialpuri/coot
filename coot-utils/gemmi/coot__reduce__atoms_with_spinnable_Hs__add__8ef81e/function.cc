#include "function.hh"

#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

gemmi::Atom *cid_to_atom_gemmi(const std::string &cid, gemmi::Structure &st) {
    gemmi::Atom *atom_p = nullptr;

    if (cid.empty()) {
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                for (auto &residue : chain.residues) {
                    if (!residue.atoms.empty()) {
                        return &residue.atoms[0];
                    }
                }
            }
        }
        return nullptr;
    }

    std::string rest = cid;
    if (rest.size() >= 2 && rest[0] == '/' && rest[1] == '/') {
        rest = rest.substr(2);
    }

    size_t pos1 = rest.find('/');
    if (pos1 == std::string::npos) {
        return nullptr;
    }

    std::string chain_id = rest.substr(0, pos1);
    std::string rest2 = rest.substr(pos1 + 1);

    size_t pos2 = rest2.find('/');
    std::string res_seq_num_str;
    std::string atom_name;

    if (pos2 == std::string::npos) {
        res_seq_num_str = rest2;
        atom_name = "";
    } else {
        res_seq_num_str = rest2.substr(0, pos2);
        atom_name = rest2.substr(pos2 + 1);
    }

    int res_seq_num = 0;
    try {
        res_seq_num = std::stoi(res_seq_num_str);
    } catch (...) {
        return nullptr;
    }

    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name != chain_id) {
                continue;
            }

            for (auto &residue : chain.residues) {
                if (residue.seqid.num.value != res_seq_num) {
                    continue;
                }

                if (atom_name.empty()) {
                    if (!residue.atoms.empty()) {
                        return &residue.atoms[0];
                    }
                } else {
                    for (auto &atom : residue.atoms) {
                        std::string trimmed = atom.name;
                        auto start = trimmed.find_first_not_of(' ');
                        if (start != std::string::npos) {
                            trimmed = trimmed.substr(start);
                        }
                        while (!trimmed.empty() && trimmed.back() == ' ') {
                            trimmed.pop_back();
                        }

                        if (trimmed == atom_name) {
                            return &atom;
                        }
                    }
                }
            }
        }
    }

    return nullptr;
}

} // namespace molecule_t
} // namespace coot