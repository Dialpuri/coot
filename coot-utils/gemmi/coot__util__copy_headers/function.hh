#pragma once

#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

bool copy_headers_gemmi(gemmi::Structure* st1, gemmi::Structure* st2, bool include_cryst);

}
}