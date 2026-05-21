#include "function.hh"

namespace coot {

void
set_has_hydrogen_atoms_state(const std::vector<gemmi::Atom>& atoms, bool& model_has_hydrogen_atoms) {
    bool found = false;
    for (const auto& atom : atoms) {
        if (is_hydrogen(atom)) {
            found = true;
            break;
        }
    }
    if (!found)
        model_has_hydrogen_atoms = false;
}

} // namespace coot