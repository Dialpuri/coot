#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/spacegroup.h>
#include <map>
#include <string>
#include <vector>
#include <utility>

namespace coot {
struct residue_spec_t;
}

namespace coot {
namespace util {
struct density_correlation_stats_info_t;
}
}

namespace coot {
namespace util {

std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
map_to_model_correlation_stats_per_residue_run_gemmi(
    const gemmi::Model &model,
    const std::string &chain_id,
    const clipper::Xmap<float> &xmap,
    unsigned int n_residues_per_blob,
    bool exclude_CON,
    float atom_mask_radius,
    float NOC_mask_radius);

} // namespace util
} // namespace coot