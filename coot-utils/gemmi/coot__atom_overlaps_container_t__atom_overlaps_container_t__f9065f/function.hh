#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {

/* ── gemmi-compatible overlap entry ── */
struct atom_overlap_t {
    int ligand_atom_index{-1};
    double overlap_volume{-1.0};
    bool is_h_bond{false};

    atom_overlap_t() = default;
};

/* ── gemmi-compatible dot color info (minimal) ── */
struct dot_color_info_t {
    double value{0.0};
};

/* ── gemmi-compatible dots container ── */
class atom_overlaps_dots_container_t {
public:
    std::vector<dot_color_info_t> dots;
    std::vector<dot_color_info_t> clashes;

    atom_overlaps_dots_container_t();

    float score() const;
};

/* ── gemmi-compatible overlaps container ── */
class atom_overlaps_container_t {
public:
    std::vector<atom_overlap_t> overlaps;

    // Ported default constructor (static factory to match naming convention)
    static atom_overlaps_container_t atom_overlaps_container_t_gemmi();

    atom_overlaps_container_t() = default;

    bool get_have_dictionary() const { return have_dictionary_; }

private:
    bool have_dictionary_{false};
};

// ── inline definitions ──

inline atom_overlaps_dots_container_t::atom_overlaps_dots_container_t()
    : dots(6), clashes() {}

inline float atom_overlaps_dots_container_t::score() const {
    return 0.0f;
}

inline atom_overlaps_container_t
atom_overlaps_container_t::atom_overlaps_container_t_gemmi() {
    return atom_overlaps_container_t{};
}

} // namespace coot
