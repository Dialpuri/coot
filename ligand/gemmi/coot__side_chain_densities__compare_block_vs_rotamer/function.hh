#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <string>
#include <utility>
#include <map>
#include <tuple>
#include <vector>

namespace coot {

struct density_box_t {
    gemmi::Residue* residue_p;
};

inline std::vector<std::string> split_string_no_blanks(const std::string& line) {
    std::vector<std::string> words;
    std::string word;
    for (char c : line) {
        if (c == ' ' || c == '\t') {
            if (!word.empty()) {
                words.push_back(word);
                word.clear();
            }
        } else {
            word += c;
        }
    }
    if (!word.empty()) {
        words.push_back(word);
    }
    return words;
}

inline int string_to_int(const std::string& s) {
    return std::stoi(s);
}

inline double string_to_double(const std::string& s) {
    return std::stod(s);
}

inline std::vector<std::string> glob_files(const std::string& dir, const std::string& pattern) {
    return std::vector<std::string>();
}

inline void get_results_addition_lock() {}
inline void release_results_addition_lock() {}

std::pair<bool, double>
compare_block_vs_rotamer_gemmi(const density_box_t& block,
                               gemmi::Residue* residue_p,
                               const std::string& rotamer_dir,
                               const clipper::Xmap<float>& xmap);

double get_log_likelihood_ratio(unsigned int grid_idx,
                                const density_box_t& block,
                                const std::string& rotamer_dir,
                                double step_size,
                                double mean,
                                double var,
                                double skew);

} // namespace coot