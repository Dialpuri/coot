#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

void find_out_of_register_errors_gemmi(gemmi::Structure* mutations_st, gemmi::Structure* ref_st);

} // namespace coot