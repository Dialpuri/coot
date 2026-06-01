#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/unitcell.hpp>
#include <vector>
#include <set>

namespace CXXCreator {

struct CXXCreatorResult {
    std::vector<const gemmi::Atom*> SelAtom;
    int nSelAtoms;
    double gridSpacing;
    double probeRadius;
};

// gemmi port of CXXCreator::CXXCreator constructor
// When context_atoms is empty (context_selHnd == 0), returns sel_atoms as-is.
// When context_atoms is non-empty, finds atoms within 35 Å of sel_atoms
// (excluding sel_atoms themselves), then intersects with context_atoms.
CXXCreatorResult CXXCreator_gemmi(gemmi::Model& model,
                                   const gemmi::UnitCell& cell,
                                   const std::vector<const gemmi::Atom*>& sel_atoms,
                                   const std::vector<const gemmi::Atom*>& context_atoms);

inline CXXCreatorResult CXXCreator_gemmi(gemmi::Model& model,
                                          const gemmi::UnitCell& cell,
                                          const std::vector<const gemmi::Atom*>& sel_atoms,
                                          const std::vector<const gemmi::Atom*>& context_atoms) {
    CXXCreatorResult result;
    result.gridSpacing = 1.2;
    result.probeRadius = 1.6;
    result.SelAtom = {};
    result.nSelAtoms = 0;

    // context_selHnd == 0 → just use the selection directly
    if (context_atoms.empty()) {
        result.SelAtom = sel_atoms;
        result.nSelAtoms = static_cast<int>(sel_atoms.size());
        return result;
    }

    // Build sets for fast lookup
    std::set<const gemmi::Atom*> sel_set(sel_atoms.begin(), sel_atoms.end());
    std::set<const gemmi::Atom*> context_set(context_atoms.begin(), context_atoms.end());

    // Find neighbors within 35 Å of sel_atoms, excluding sel_atoms themselves
    // (matching MMDB SelectNeighbours behaviour)
    gemmi::NeighborSearch ns(model, cell, 35.0);
    ns.populate(false);

    std::set<const gemmi::Atom*> neighbor_set;
    for (const auto* atom : sel_atoms) {
        auto hits = ns.find_atoms(atom->pos, '\0', 0.0, 35.0);
        for (auto* mark : hits) {
            gemmi::CRA cra = mark->to_cra(model);
            if (cra.atom && !sel_set.count(cra.atom)) {
                neighbor_set.insert(cra.atom);
            }
        }
    }

    // AND with context selection
    for (const auto* atom : neighbor_set) {
        if (context_set.count(atom)) {
            result.SelAtom.push_back(atom);
        }
    }

    result.nSelAtoms = static_cast<int>(result.SelAtom.size());
    return result;
}

} // namespace CXXCreator