#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace ideal_rna {

void add_o2_prime_gemmi(gemmi::Residue* res);

}
}