#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <vector>

namespace coot {
namespace util {

void sfcalc_genmap_gemmi(
    gemmi::Model *mol,
    const clipper::HKL_data<clipper::data32::F_sigF> &fobs_in,
    const clipper::HKL_data<clipper::data32::Flag> &free,
    clipper::Xmap<float> *xmap_p);

} // namespace util
} // namespace coot