#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>
#include "coot/coot-utils/coot-shelx.hh"

namespace coot {

struct shelx_read_file_info_t_gemmi {
    int status;
    int udd_afix_handle;
    std::shared_ptr<gemmi::Structure> structure;
    int is_protein_flag;
    
    shelx_read_file_info_t_gemmi() : status(0), udd_afix_handle(-1), structure(nullptr), is_protein_flag(0) {}
    shelx_read_file_info_t_gemmi(int st, int udd, std::shared_ptr<gemmi::Structure> str)
        : status(st), udd_afix_handle(udd), structure(str), is_protein_flag(0) {}
    shelx_read_file_info_t_gemmi(int st, int udd, std::shared_ptr<gemmi::Structure> str, int protein)
        : status(st), udd_afix_handle(udd), structure(str), is_protein_flag(protein) {}
};

inline shelx_read_file_info_t_gemmi
read_file_gemmi(ShelxIns& s, const std::string& filename) {
    shelx_read_file_info_t_gemmi result;
    result.status = 0;
    result.udd_afix_handle = -1;
    result.structure = nullptr;
    result.is_protein_flag = 0;
    
    std::ifstream f(filename.c_str());
    if (!f) {
        return result;
    }
    
    result.status = 1;
    result.structure = std::make_shared<gemmi::Structure>();
    result.structure->name = "ShelxIns";
    
    gemmi::Model model;
    gemmi::Chain chain("A");
    
    std::vector<gemmi::Atom> atom_vector;
    std::string current_res_name = "";
    int current_res_no = 0;
    int resi_count = 0;
    int udd_afix_handle = 1;
    
    std::string line;
    while (std::getline(f, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) {
            line.pop_back();
        }
        
        if (line.empty()) continue;
        
        std::vector<std::string> words;
        std::string current_word;
        bool in_quotes = false;
        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            if (c == '"') {
                in_quotes = !in_quotes;
            } else if (!in_quotes && (c == ' ' || c == '\t')) {
                if (!current_word.empty()) {
                    words.push_back(current_word);
                    current_word.clear();
                }
            } else {
                current_word += c;
            }
        }
        if (!current_word.empty()) {
            words.push_back(current_word);
        }
        
        if (words.empty()) continue;
        
        std::string card_word_0 = words[0];
        for (char& c : card_word_0) c = std::toupper(c);
        
        if (card_word_0 == "RESI") {
            if (!atom_vector.empty()) {
                gemmi::Residue res;
                res.name = current_res_name;
                res.seqid.num.value = current_res_no;
                res.atoms = atom_vector;
                chain.residues.push_back(res);
                atom_vector.clear();
            }
            resi_count++;
            
            if (words.size() > 1) {
                std::string res_no_string = words[1];
                size_t p = res_no_string.find(':');
                if (p != std::string::npos) {
                    std::string res_chain_id = res_no_string.substr(0, p);
                    std::string r_string = res_no_string.substr(p + 1);
                    if (res_chain_id != "A") {
                        chain.name = res_chain_id;
                    }
                    try {
                        current_res_no = std::stoi(r_string);
                    } catch (...) {}
                } else {
                    try {
                        current_res_no = std::stoi(res_no_string);
                    } catch (...) {}
                }
            }
            
            if (words.size() > 2) {
                current_res_name = words[2];
            } else {
                current_res_name = "";
            }
        } else if (card_word_0 == "END" || card_word_0 == "HKLF") {
            if (!atom_vector.empty()) {
                gemmi::Residue res;
                res.name = current_res_name;
                res.seqid.num.value = current_res_no;
                res.atoms = atom_vector;
                chain.residues.push_back(res);
                atom_vector.clear();
            }
            if (card_word_0 == "HKLF") break;
        } else if (card_word_0 == "CELL") {
            if (words.size() >= 7) {
                double a = 10.0, b = 10.0, c = 10.0;
                double alpha = 90.0, beta = 90.0, gamma = 90.0;
                try {
                    a = std::stod(words[2]);
                    b = std::stod(words[3]);
                    c = std::stod(words[4]);
                    alpha = std::stod(words[5]);
                    beta = std::stod(words[6]);
                    gamma = std::stod(words[7]);
                } catch (...) {}
                result.structure->cell.set(a, b, c, alpha, beta, gamma);
            }
        } else if (words.size() > 4) {
            try {
                gemmi::Atom atom;
                atom.name = words[0];
                int coord_start = 2;
                if (words.size() >= coord_start + 3) {
                    atom.pos.x = std::stod(words[coord_start]);
                    atom.pos.y = std::stod(words[coord_start + 1]);
                    atom.pos.z = std::stod(words[coord_start + 2]);
                }
                if (words.size() > coord_start + 3) {
                    atom.occ = std::stod(words[coord_start + 3]);
                }
                if (words.size() > coord_start + 4) {
                    atom.b_iso = std::stod(words[coord_start + 4]);
                }
                if (words.size() > coord_start + 5) {
                    atom.element = gemmi::Element(words[coord_start + 5][0]);
                }
                atom_vector.push_back(atom);
            } catch (...) {}
        }
    }
    
    if (!atom_vector.empty()) {
        gemmi::Residue res;
        res.name = current_res_name;
        res.seqid.num.value = current_res_no;
        res.atoms = atom_vector;
        chain.residues.push_back(res);
    }
    
    model.chains.push_back(chain);
    result.structure->models.push_back(model);
    result.udd_afix_handle = udd_afix_handle;
    
    if (resi_count > 10) {
        result.is_protein_flag = 1;
    }
    
    return result;
}

} // namespace coot