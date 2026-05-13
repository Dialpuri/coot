#pragma once

#include <map>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot { namespace util {

std::map<std::string, clipper::RTop_orth>
get_ori_to_this_res_gemmi(const gemmi::Residue& residue);

}}