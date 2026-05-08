#include "function.hh"

namespace coot {
namespace molecule_t {

gemmi::Atom *cid_to_atom_gemmi(const std::string &cid, gemmi::Structure &st) {
    gemmi::Atom *atom_p = nullptr;

    if (cid.empty()) {
        // Empty CID returns first atom in the structure
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

    // Parse CID: "//chain/residue_num/atom_name"
    // Skip the leading "//"
    std::string rest = cid;
    if (rest.size() >= 2 && rest[0] == '/' && rest[1] == '/') {
        rest = rest.substr(2);
    }

    // Split by '/'
    size_t pos1 = rest.find('/');
    if (pos1 == std::string::npos) {
        // No slash found - invalid format
        return nullptr;
    }

    std::string chain_id = rest.substr(0, pos1);
    std::string rest2 = rest.substr(pos1 + 1);

    size_t pos2 = rest2.find('/');
    std::string res_seq_num_str;
    std::string atom_name;

    if (pos2 == std::string::npos) {
        // Only chain/residue, no atom name
        res_seq_num_str = rest2;
        atom_name = "";
    } else {
        res_seq_num_str = rest2.substr(0, pos2);
        atom_name = rest2.substr(pos2 + 1);
    }

    // Convert residue sequence number to int
    int res_seq_num = 0;
    try {
        res_seq_num = std::stoi(res_seq_num_str);
    } catch (...) {
        return nullptr;
    }

    // Search for the atom
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
                    // No atom name specified, return first atom in residue
                    if (!residue.atoms.empty()) {
                        return &residue.atoms[0];
                    }
                } else {
                    // Look for specific atom name
                    for (auto &atom : residue.atoms) {
                        // Trim the atom name (PDB format has padded names)
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