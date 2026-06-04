#pragma once

#include <clipper/core/xmap.h>
#include <clipper/core/cell.h>
#include <clipper/core/symop.h>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <gemmi/resinfo.hpp>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <cmath>
#include <iostream>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

std::vector<std::pair<std::string, clipper::Xmap<float>>>
partition_map_by_chain_gemmi(const clipper::Xmap<float>& xmap,
                             const gemmi::Model* model_p,
                             std::string* state_string_p);

} // namespace util
} // namespace coot