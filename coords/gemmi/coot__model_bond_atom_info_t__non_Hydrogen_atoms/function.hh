#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

class model_bond_atom_info_t_gemmi {
public:
    void add_bond_atom(const gemmi::Atom& atom) {
        atoms_.push_back(&atom);
        if (atom.element.is_hydrogen()) {
            hydrogen_atoms_.push_back(&atom);
        } else {
            non_hydrogen_atoms_.push_back(&atom);
        }
    }

    int n_H() const { return static_cast<int>(hydrogen_atoms_.size()); }

    int n_non_H() const { return static_cast<int>(non_hydrogen_atoms_.size()); }

    std::vector<const gemmi::Atom*> Hydrogen_atoms() const {
        return hydrogen_atoms_;
    }

    std::vector<const gemmi::Atom*> non_Hydrogen_atoms_gemmi() const {
        return non_hydrogen_atoms_;
    }

private:
    std::vector<const gemmi::Atom*> atoms_;
    std::vector<const gemmi::Atom*> hydrogen_atoms_;
    std::vector<const gemmi::Atom*> non_hydrogen_atoms_;
};

} // namespace coot