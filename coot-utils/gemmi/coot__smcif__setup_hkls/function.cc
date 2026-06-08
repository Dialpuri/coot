#include "function.hh"

#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>

namespace coot {
namespace smcif {

std::vector<clipper::HKL> setup_hkls_gemmi(const std::string& file_name) {
    std::vector<clipper::HKL> hkls;
    
    gemmi::cif::Document doc;
    try {
        doc = gemmi::read_cif_gz(file_name);
    } catch (...) {
        std::cout << "WARNING:: Error reading small-molecule cif \"" << file_name
                  << "\"" << std::endl;
        return hkls;
    }
    
    if (doc.blocks.empty()) {
        std::cout << "WARNING:: Error reading small-molecule cif \"" << file_name
                  << "\"" << std::endl;
        return hkls;
    }
    
    gemmi::cif::Block& block = doc.blocks[0];
    
    // Try standard refln loop first
    std::string h_tag = "_refln_index_h";
    std::string k_tag = "_refln_index_k";
    std::string l_tag = "_refln_index_l";
    
    gemmi::cif::Column h_col = block.find_values(h_tag);
    if (!h_col) {
        // Try pd_refln loop as fallback
        h_tag = "_pd_refln_index_h";
        k_tag = "_pd_refln_index_k";
        l_tag = "_pd_refln_index_l";
        h_col = block.find_values(h_tag);
    }
    
    if (h_col) {
        gemmi::cif::Column k_col = block.find_values(k_tag);
        gemmi::cif::Column l_col = block.find_values(l_tag);
        
        if (!k_col || !l_col) return hkls;
        
        int length = h_col.length();
        
        for (int i = 0; i < length; ++i) {
            std::string h_str = h_col[i];
            std::string k_str = k_col[i];
            std::string l_str = l_col[i];
            
            if (h_str.empty() || k_str.empty() || l_str.empty()) continue;
            if (gemmi::cif::is_null(h_str) || gemmi::cif::is_null(k_str) || gemmi::cif::is_null(l_str)) continue;
            
            int h, k, l;
            try {
                h = std::stoi(h_str);
                k = std::stoi(k_str);
                l = std::stoi(l_str);
            } catch (...) {
                continue;
            }
            
            hkls.push_back(clipper::HKL(h, k, l));
        }
    }
    
    return hkls;
}

}  // namespace smcif
}  // namespace coot