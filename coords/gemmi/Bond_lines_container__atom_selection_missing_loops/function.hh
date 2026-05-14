#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

struct atom_selection_container_t {
    gemmi::Structure* structure;
    gemmi::Model* model;
    int n_selected_atoms;
    int UDDAtomIndexHandle;
    int mol;
};

class Bond_lines_container {
public:
    void atom_selection_missing_loops_gemmi(const atom_selection_container_t &asc,
                                            int udd_atom_index_handle,
                                            int udd_fixed_during_refinement_handle);
};

} // namespace coot