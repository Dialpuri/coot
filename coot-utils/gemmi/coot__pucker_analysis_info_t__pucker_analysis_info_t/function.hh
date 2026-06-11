#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <stdexcept>
#include <cmath>

namespace coot {

// ---- lsq_plane_info_t (inline re-implementation using gemmi::Vec3) ----

class lsq_plane_info_t {
public:
    gemmi::Vec3 normal_;
    double d_;

    lsq_plane_info_t() : normal_(0, 0, 0), d_(0) {}

    explicit lsq_plane_info_t(const std::vector<gemmi::Vec3>& v);

    gemmi::Vec3 centre() const;

    gemmi::Vec3 projected_point(gemmi::Vec3 p) const;
};

// ---- pucker_analysis_markup_info_t ----

struct pucker_analysis_markup_info_t {
    gemmi::Vec3 base_ring_centre;
    gemmi::Vec3 base_ring_normal;
    gemmi::Vec3 phosphorus_position;
    gemmi::Vec3 projected_point;
};

// ---- pucker_analysis_info_t ----

class pucker_analysis_info_t {
public:
    enum PUCKERED_ATOM_T {
        NONE = 0,
        C1_PRIME = 1,
        C2_PRIME = 2,
        C3_PRIME = 3,
        C4_PRIME = 4,
        O4_PRIME = 5
    };

    std::vector<gemmi::Vec3> base_atoms_coords;
    std::vector<gemmi::Vec3> ribose_atoms_coords;
    double out_of_plane_distance;
    double plane_distortion;
    gemmi::Atom* C1_prime;
    gemmi::Atom* N1_or_9;
    std::string altconf;
    pucker_analysis_markup_info_t markup_info;
    PUCKERED_ATOM_T puckered_atom_;

    pucker_analysis_info_t()
        : out_of_plane_distance(0), plane_distortion(0),
          C1_prime(nullptr), N1_or_9(nullptr), puckered_atom_(NONE) {}

    // The gemmi ported constructor equivalent
    void pucker_analysis_info_t_gemmi(gemmi::Residue const& res, std::string altconf_in);
};

// ---- lsq_plane_deviation (free function) ----

std::pair<double, double> lsq_plane_deviation(
    const std::vector<gemmi::Vec3>& plane_atom_coords,
    gemmi::Vec3 pt);

} // namespace coot
