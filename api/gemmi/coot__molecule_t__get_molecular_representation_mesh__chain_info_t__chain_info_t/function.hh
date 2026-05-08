#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {
namespace get_molecular_representation_mesh {

struct chain_info_t {
    gemmi::Chain* chain_p;
    int resno_min;
    int resno_max;
};

inline chain_info_t chain_info_t_gemmi(gemmi::Chain* c, int min, int max) {
    return chain_info_t{c, min, max};
}

} // namespace get_molecular_representation_mesh
} // namespace molecule_t
} // namespace coot