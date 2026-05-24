#include "function.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <clipper/core/coords.h>

#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

// Helper function to convert string to int
static int string_to_int(const std::string& str) {
    size_t pos;
    int val = std::stoi(str, &pos);
    if (pos != str.size()) {
        throw std::runtime_error("Invalid integer string");
    }
    return val;
}

// Helper function to convert string to float
static float string_to_float(const std::string& str) {
    size_t pos;
    float val = std::stof(str, &pos);
    if (pos != str.size()) {
        throw std::runtime_error("Invalid float string");
    }
    return val;
}

// Helper function to trim whitespace
static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<coot::cablam_markup_t>
coot::make_cablam_markups_gemmi(gemmi::Structure &st, const std::string &cablam_output_file_name) {
    std::vector<cablam_markup_t> v;

    std::vector<std::pair<residue_spec_t, double> > scored_baddie_specs;
    
    std::ifstream f(cablam_output_file_name.c_str());
    if (f) {
        std::string line;
        while (std::getline(f, line)) {
            // Parse by position, not field!
            if (line.length() == 90) {
                std::string chain_id = line.substr(0,2);
                if (chain_id[0] == ' ') chain_id = line.substr(1,1);
                chain_id = trim(chain_id);
                
                std::string resno_string = line.substr(2,4);
                try {
                    int res_no = string_to_int(resno_string);
                    std::string residue_type = line.substr(8,3);
                    residue_type = trim(residue_type);
                    
                    std::string cablam_string = line.substr(13,6);
                    cablam_string = trim(cablam_string);
                    
                    std::string cablam_type_string = line.substr(20,7);
                    cablam_type_string = trim(cablam_type_string);
                    
                    if (cablam_string == "CaBLAM") {
                        // either Disfavoured or an Outlier
                        if (cablam_type_string == "Outlier") {
                            std::string level_string = line.substr(34,6);
                            double level = string_to_float(level_string);

                            residue_spec_t res_spec(chain_id, res_no, "");
                            std::pair<residue_spec_t, double> p(res_spec, level);
                            scored_baddie_specs.push_back(p);
                        }
                    }
                }
                catch (const std::runtime_error &rte) {
                    std::cout << "Error parsing: " << rte.what() << std::endl;
                }
            }
        }
        f.close();
    } else {
        std::cout << "WARNING:: file not found " << cablam_output_file_name << std::endl;
    }
    
    v = make_cablam_markups_gemmi(scored_baddie_specs, st);
    return v;
}

std::vector<coot::cablam_markup_t>
coot::make_cablam_markups_gemmi(const std::vector<std::pair<residue_spec_t, double> > &residues,
                                gemmi::Structure &st) {
    std::vector<cablam_markup_t> v;
    
    // Find residues in the structure that match the specs
    for (const auto& spec_pair : residues) {
        const residue_spec_t& spec = spec_pair.first;
        double score = spec_pair.second;
        
        // Search for matching residue in structure
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& residue : chain.residues) {
                    // Check if this residue matches the spec
                    bool chain_match = (chain.name == spec.chain_id);
                    bool resno_match = (residue.seqid.num.value == spec.res_no);
                    
                    // For insertion code, spec uses "" for none, gemmi uses ' '
                    char icode = residue.seqid.icode;
                    bool icode_match = (spec.ins_code.empty() || spec.ins_code == " " || spec.ins_code == std::string(1, icode));
                    
                    if (chain_match && resno_match && icode_match) {
                        cablam_markup_t markup;
                        markup.score = score;
                        markup.cra = {&chain, &residue, nullptr};
                        v.push_back(markup);
                        break; // Found match, move to next spec
                    }
                }
            }
        }
    }
    
    return v;
}