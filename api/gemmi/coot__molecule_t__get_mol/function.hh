#pragma once
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct molecule_t {
    gemmi::Structure structure;
    
    // Original: mmdb::Manager *get_mol() { return atom_sel.mol; }
    // Ported: return pointer to underlying gemmi Structure
    gemmi::Structure* get_mol_gemmi() { return &structure; }
};

} // namespace coot