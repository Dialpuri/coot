#pragma once
#include <gemmi/model.hpp>
#include "coords/Cartesian.hh"

struct molecule_extents_t {
    coot::Cartesian left;
    coot::Cartesian right;
    coot::Cartesian front;
    coot::Cartesian back;
    coot::Cartesian bottom;
    coot::Cartesian top;
    coot::Cartesian centre;
    float expansion_size_;

    molecule_extents_t() : expansion_size_(0.0f) {}

    /// gemmi port of molecule_extents_t constructor.
    /// Takes a gemmi::Model instead of atom_selection_container_t.
    static molecule_extents_t molecule_extents_t_gemmi(
        const gemmi::Model& model, float expansion_size);
};