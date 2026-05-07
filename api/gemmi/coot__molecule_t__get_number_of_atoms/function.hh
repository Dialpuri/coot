#pragma once
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {
    unsigned int get_number_of_atoms_gemmi(const gemmi::Structure& st);
}
}