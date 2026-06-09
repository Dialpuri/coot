#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

struct density_box_t {
    float *density_box;
    gemmi::Residue *residue_p;
    int n_steps;
    double mean;
    double var;
    bool is_weird;

    void init() {
        density_box = nullptr;
        residue_p = nullptr;
        n_steps = 0;
        mean = 0.;
        var = -1.;
        is_weird = false;
    }

    density_box_t() { init(); }
    density_box_t(float *density_box, gemmi::Residue *residue_p, int n_steps);

    // gemmi ported factory methods
    static density_box_t density_box_t_gemmi();
    static density_box_t density_box_t_gemmi(float *density_box, gemmi::Residue *residue_p, int n_steps);

    density_box_t(const density_box_t& other) {
        density_box = other.density_box;
        residue_p = other.residue_p;
        n_steps = other.n_steps;
        mean = other.mean;
        var = other.var;
        is_weird = other.is_weird;
    }

    density_box_t &operator=(const density_box_t& other) {
        density_box = other.density_box;
        residue_p = other.residue_p;
        n_steps = other.n_steps;
        mean = other.mean;
        var = other.var;
        is_weird = other.is_weird;
        return *this;
    }

    bool empty() const { return density_box == nullptr; }
    int nnn() const {
        int n = 2 * n_steps + 1;
        return n * n * n;
    }
    float &operator[](int i) {
        return density_box[i];
    }
    const float &operator[](int i) const {
        return density_box[i];
    }
    float *box() { return density_box; }
};

// Inline definitions for short factory methods
inline density_box_t density_box_t::density_box_t_gemmi() {
    density_box_t box;
    box.init();
    return box;
}

inline density_box_t density_box_t::density_box_t_gemmi(float *density_box, gemmi::Residue *residue_p, int n_steps) {
    density_box_t box;
    box.init();
    box.density_box = density_box;
    box.residue_p = residue_p;
    box.n_steps = n_steps;
    return box;
}

} // namespace coot
