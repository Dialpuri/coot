#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {

class my_atom_colour_map_t {
public:
    std::vector<std::string> atom_colour_map;

    my_atom_colour_map_t() {
        for (int i = 0; i < 50; i++) {
            atom_colour_map.push_back("---");
        }
    }

    int index_for_chain(const std::string &chain_id) {
        int n = static_cast<int>(atom_colour_map.size());
        int idx = -1;
        for (int i = 0; i < n; i++) {
            if (atom_colour_map[i] == chain_id) {
                idx = i;
                break;
            }
        }
        if (idx == -1) {
            atom_colour_map.push_back(chain_id);
            idx = static_cast<int>(atom_colour_map.size()) - 1;
        }
        return idx;
    }

    void fill_chain_id_map_gemmi(const gemmi::Model &model) {
        for (const gemmi::Chain &chain : model.chains) {
            if (!chain.residues.empty()) {
                std::string chain_id = chain.name;
                unsigned int idx = static_cast<unsigned int>(index_for_chain(chain_id));
            }
        }
    }
};

} // namespace coot