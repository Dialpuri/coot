#include "function.hh"

namespace coot {

double get_log_likelihood_ratio(unsigned int grid_idx,
                                const density_box_t& block,
                                const std::string& rotamer_dir,
                                double step_size,
                                double mean,
                                double var,
                                double skew) {
    return 0.0;
}

std::pair<bool, double>
compare_block_vs_rotamer_gemmi(const density_box_t& block,
                               gemmi::Residue* residue_p,
                               const std::string& rotamer_dir,
                               const clipper::Xmap<float>& xmap) {
    bool success = false;
    double sum_log_likelihood = 0.0;
    
    get_results_addition_lock();
    
    std::map<std::string, std::map<unsigned int, std::tuple<double, double, double>>> rotamer_dir_grid_stats_map_cache;
    
    auto it = rotamer_dir_grid_stats_map_cache.find(rotamer_dir);
    
    if (it == rotamer_dir_grid_stats_map_cache.end()) {
        std::map<unsigned int, std::tuple<double, double, double>> stats_map;
        rotamer_dir_grid_stats_map_cache[rotamer_dir] = stats_map;
        it = rotamer_dir_grid_stats_map_cache.find(rotamer_dir);
    }
    
    if (it != rotamer_dir_grid_stats_map_cache.end()) {
        const std::map<unsigned int, std::tuple<double, double, double>>& stats_map = it->second;
        if (!stats_map.empty()) {
            success = true;
            for (const auto& pair : stats_map) {
                const unsigned int& grid_idx = pair.first;
                const std::tuple<double, double, double>& m_v_s = pair.second;
                const double& mean = std::get<0>(m_v_s);
                const double& var = std::get<1>(m_v_s);
                const double& skew = std::get<2>(m_v_s);
                
                double llr = get_log_likelihood_ratio(grid_idx, block, rotamer_dir, 0.0, mean, var, skew);
                sum_log_likelihood += llr;
            }
        }
    }
    
    release_results_addition_lock();
    
    return std::make_pair(success, sum_log_likelihood);
}

} // namespace coot