#pragma once

#include <clipper/clipper.h>
#include <gemmi/model.hpp>

namespace coot {

void shift_field_xyz_refinement_gemmi(
    const clipper::HKL_data<clipper::datatypes::F_sigF<float>> &fo0,
    const clipper::HKL_data<clipper::data32::Flag> &free,
    gemmi::Structure &mol,
    float resolution) {
}

}  // namespace coot