#pragma once
#include <vector>
#include <iostream>
#include <gemmi/model.hpp>

namespace coot {
namespace util {
namespace map_to_model_correlation_stats_per_residue_run {

class residue_run_t {
public:
    unsigned int idx_mid;
    unsigned int n_residues_per_blob;
    std::vector<gemmi::Residue> residues;
    
    residue_run_t() : idx_mid(0), n_residues_per_blob(0) {}
    explicit residue_run_t(unsigned int n_residues_per_blob_in) 
        : idx_mid(n_residues_per_blob_in/2), n_residues_per_blob(n_residues_per_blob_in) {}
    explicit residue_run_t(const std::vector<gemmi::Residue> &rr_in) 
        : residues(rr_in), n_residues_per_blob(rr_in.size()), idx_mid(rr_in.size()/2) {}
    
    const gemmi::Residue* residue_mid() const {
        if (residues.size() >= n_residues_per_blob) {
            return &residues[idx_mid];
        } else {
            std::cout << "ERROR:: indexing residues_mid " << residues.size() << std::endl;
            return nullptr;
        }
    }
    
    void add_residue(const gemmi::Residue& r) { residues.push_back(r); }
    void add(const std::vector<gemmi::Residue> &rv) {
        residues.insert(residues.begin(), rv.cbegin(), rv.cend());
    }
};

} // namespace map_to_model_correlation_stats_per_residue_run
} // namespace util
} // namespace coot