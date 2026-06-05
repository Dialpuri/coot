#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

struct elastic_network_item_t_gemmi {
    const gemmi::Atom* at_1;
    const gemmi::Atom* at_2;
    float spring_constant;

    elastic_network_item_t_gemmi()
        : at_1(nullptr), at_2(nullptr), spring_constant(0.0f) {}

    elastic_network_item_t_gemmi(const gemmi::Atom* a1,
                                  const gemmi::Atom* a2,
                                  float k)
        : at_1(a1), at_2(a2), spring_constant(k) {}
};

struct elastic_network_model_t_gemmi {
    std::vector<elastic_network_item_t_gemmi> d;

    elastic_network_model_t_gemmi() = default;

    elastic_network_model_t_gemmi(
        const std::vector<const gemmi::Atom*>& selected_atoms,
        double min_dist,
        double max_dist,
        unsigned int max_n_distances);
};

} // namespace coot