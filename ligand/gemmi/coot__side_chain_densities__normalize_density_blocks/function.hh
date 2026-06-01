#pragma once

#include <gemmi/model.hpp>
#include <map>

namespace coot {

struct density_box_t_gemmi {
    float* density_box;
    gemmi::Residue* residue_p;
    int n_steps;

    density_box_t_gemmi() : density_box(nullptr), residue_p(nullptr), n_steps(0) {}
    density_box_t_gemmi(float* box, gemmi::Residue* res, int steps)
        : density_box(box), residue_p(res), n_steps(steps) {}

    // Copy constructor - shallow copy (same as original density_box_t)
    density_box_t_gemmi(const density_box_t_gemmi& other) = default;
    density_box_t_gemmi& operator=(const density_box_t_gemmi& other) = default;

    bool empty() const { return n_steps == 0; }

    int nnn() const {
        int n = 2 * n_steps + 1;
        return n * n * n;
    }

    float operator[](unsigned int idx) const {
        return density_box[idx];
    }

    void scale_by(float scale_factor) {
        if (n_steps > 0) {
            int n = 2 * n_steps + 1;
            int n3 = n * n * n;
            for (int i = 0; i < n3; i++)
                if (density_box[i] > -1000)
                    density_box[i] *= scale_factor;
        }
    }
};

namespace side_chain_densities {

inline void normalize_density_blocks_gemmi(
    std::map<gemmi::Residue*, density_box_t_gemmi>& density_block_map_cache,
    double mn_scale_for_normalized_density) {

    unsigned int n_grid_pts = 0;
    double sum = 0;

    for (auto it = density_block_map_cache.begin(); it != density_block_map_cache.end(); ++it) {
        const density_box_t_gemmi& block = it->second;
        if (!block.empty()) {
            int nnn = block.nnn();
            for (int i = 0; i < nnn; i++) {
                if (block[i] > 0.0) {
                    sum += block[i];
                    n_grid_pts++;
                }
            }
        }
    }

    if (n_grid_pts > 0) {
        double av = sum / static_cast<double>(n_grid_pts);
        double sc = mn_scale_for_normalized_density / av;
        for (auto it = density_block_map_cache.begin(); it != density_block_map_cache.end(); ++it) {
            density_box_t_gemmi& block = it->second;
            block.scale_by(static_cast<float>(sc));
        }
    }
}

} // namespace side_chain_densities

} // namespace coot