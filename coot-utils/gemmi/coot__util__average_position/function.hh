#pragma once
#include <vector>
#include <clipper/clipper.h>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

inline clipper::Coord_orth average_position_gemmi(const std::vector<clipper::Coord_orth>& pts) {
    if (pts.size() > 0) {
        double xsum = 0.0;
        double ysum = 0.0;
        double zsum = 0.0;
        for (size_t i = 0; i < pts.size(); ++i) {
            xsum += pts[i].x();
            ysum += pts[i].y();
            zsum += pts[i].z();
        }
        double denom = 1.0 / static_cast<double>(pts.size());
        return clipper::Coord_orth(denom * xsum, denom * ysum, denom * zsum);
    } else {
        return clipper::Coord_orth(0, 0, 0);
    }
}

} // namespace util
} // namespace coot