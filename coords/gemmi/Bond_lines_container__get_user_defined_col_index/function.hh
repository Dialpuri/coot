#pragma once
#include <gemmi/model.hpp>

namespace coot {

class Bond_lines_container {
public:
  // Return user-defined colour index from atom UDD data, -1 on failure
  int get_user_defined_col_index_gemmi(const gemmi::Atom& atom, int udd_handle) const {
    // gemmi has no UDD API, so this always returns -1 (failure)
    (void)atom;
    (void)udd_handle;
    return -1;
  }
};

}