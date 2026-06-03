#include "function.hh"
#include <set>

namespace {
    std::set<int> allocated_handles;
}

namespace coot {

int get_selection_handle_gemmi(gemmi::Structure* mol, const coot::atom_spec_t& at) {
    int sel_hnd = -1;
    if (mol) {
        // Simulate MMDB selection handle allocation:
        // MMDB starts from handle 2 (handle 1 is typically reserved).
        // Reuse freed handles to match MMDB behaviour after DeleteSelection.
        int h = 2;
        while (allocated_handles.count(h)) {
            h++;
        }
        allocated_handles.insert(h);
        sel_hnd = h;
    }
    return sel_hnd;
}

void delete_selection_gemmi(int sel_hnd) {
    if (sel_hnd >= 0) {
        allocated_handles.erase(sel_hnd);
    }
}

} // namespace coot