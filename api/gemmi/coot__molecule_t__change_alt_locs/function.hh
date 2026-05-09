#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot { namespace molecule_t {

// Split by delimiter, preserving leading empty parts (for CID "//A/1" → ["", "", "A", "1"])
inline std::vector<std::string> split_string(const std::string& s, const std::string& delimiter) {
    std::vector<std::string> parts;
    size_t start = 0;
    size_t end = s.find(delimiter);
    while (end != std::string::npos) {
        std::string part = s.substr(start, end - start);
        parts.push_back(part);
        start = end + delimiter.length();
        end = s.find(delimiter, start);
    }
    std::string last = s.substr(start);
    parts.push_back(last);
    return parts;
}

inline bool is_main_chain_p(const gemmi::Atom& at) {
    const std::string& name = at.name;
    if (name == " N  " || name == " C  " || name == " CA " ||
        name == " H  " || name == " HA " || name == " OXT" ||
        name == " O  ") {
        return true;
    }
    if (name == " CB ") return true;
    if (name == " HA2" || name == " HA3") return true;
    return false;
}

inline int swap_alt_locs_all(gemmi::Residue& res) {
    int status = 0;
    for (auto& atom : res.atoms) {
        if (atom.altloc == 'A') { atom.altloc = 'B'; status = 1; }
        else if (atom.altloc == 'B') { atom.altloc = 'A'; status = 1; }
    }
    return status;
}

inline int swap_alt_locs_sidechain(gemmi::Residue& res) {
    int status = 0;
    for (auto& atom : res.atoms) {
        if (!is_main_chain_p(atom)) {
            if (atom.altloc == 'A') { atom.altloc = 'B'; status = 1; }
            else if (atom.altloc == 'B') { atom.altloc = 'A'; status = 1; }
        }
    }
    return status;
}

inline int swap_alt_locs_mainchain(gemmi::Residue& res) {
    int status = 0;
    for (auto& atom : res.atoms) {
        if (is_main_chain_p(atom)) {
            if (atom.altloc == 'A') { atom.altloc = 'B'; status = 1; }
            else if (atom.altloc == 'B') { atom.altloc = 'A'; status = 1; }
        }
    }
    return status;
}

inline int swap_alt_locs_atomlist(gemmi::Residue& res, const std::string& change_mode) {
    int status = 0;
    auto names = split_string(change_mode, ",");
    for (const auto& name_raw : names) {
        std::string name = name_raw;
        if (name.length() == 3) name = " " + name;
        else if (name.length() == 2) name = " " + name + " ";
        else if (name.length() == 1) name = " " + name + "  ";
        for (auto& atom : res.atoms) {
            if (atom.name == name) {
                if (atom.altloc == 'A') { atom.altloc = 'B'; status = 1; }
                else if (atom.altloc == 'B') { atom.altloc = 'A'; status = 1; }
            }
        }
    }
    return status;
}

inline int change_alt_locs_gemmi(gemmi::Structure& st,
                                 const std::string& cid,
                                 const std::string& change_mode) {
    auto parts = split_string(cid, "/");
    if (parts.size() < 4) return -1;
    std::string chain_name = parts[2];
    int seqnum = std::stoi(parts[3]);

    std::vector<gemmi::Residue*> matching;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == chain_name) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == seqnum) {
                        matching.push_back(&res);
                    }
                }
            }
        }
    }

    if (matching.empty()) return -1;

    int status = 0;
    for (auto* res : matching) {
        if (change_mode == "residue") status |= swap_alt_locs_all(*res);
        else if (change_mode == "side-chain") status |= swap_alt_locs_sidechain(*res);
        else if (change_mode == "main-chain") status |= swap_alt_locs_mainchain(*res);
        else status |= swap_alt_locs_atomlist(*res, change_mode);
    }
    return status;
}

}}