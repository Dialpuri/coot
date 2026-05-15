#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

struct restraints_container_t {
    std::vector<gemmi::Connection> connections;

    void fill_links_gemmi(const gemmi::Structure* st);
};

inline void restraints_container_t::fill_links_gemmi(const gemmi::Structure* st) {
    connections.clear();
    if (st) {
        for (const auto& conn : st->connections) {
            connections.push_back(conn);
        }
    }
}

} // namespace coot