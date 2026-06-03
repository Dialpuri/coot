#pragma once

#include <vector>
#include <map>
#include <string>
#include <gemmi/model.hpp>
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"

namespace coot {

class atom_tree_t {
public:
    std::vector<atom_vertex> atom_vertex_vec;
    std::map<std::string, map_index_t> name_to_index;

    atom_tree_t(const std::vector<std::vector<int> > &contact_indices,
                int base_atom_index,
                gemmi::Residue *res,
                const std::string &altconf)
        : residue_g(res) {}

    map_index_t get_index_gemmi(gemmi::Atom *atom) const {
        map_index_t idx;
        if (residue_g) {
            for (int iat = 0; iat < static_cast<int>(residue_g->atoms.size()); iat++) {
                if (&(residue_g->atoms[iat]) == atom) {
                    idx = map_index_t(iat);
                    break;
                }
            }
        }
        return idx;
    }

protected:
    gemmi::Residue *residue_g;
};

} // namespace coot