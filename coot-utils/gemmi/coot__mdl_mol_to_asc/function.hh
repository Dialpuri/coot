#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "lidia-core/lbg-molfile.hh"

namespace coot {

struct atom_selection_container_t_gemmi {
    gemmi::Structure* mol;
    int n_selected_atoms;
    int read_success;

    atom_selection_container_t_gemmi()
        : mol(nullptr), n_selected_atoms(0), read_success(0) {}
};

inline atom_selection_container_t_gemmi mdl_mol_to_asc_gemmi(
    const lig_build::molfile_molecule_t &m)
{
    atom_selection_container_t_gemmi asc;

    if (m.atoms.empty()) {
        return asc;
    }

    gemmi::Structure* structure = new gemmi::Structure();
    gemmi::Model model;
    structure->models.push_back(std::move(model));
    gemmi::Model& mdl = structure->models.back();

    gemmi::Chain chain;
    chain.name = "A";
    mdl.chains.push_back(std::move(chain));
    gemmi::Chain& ch = mdl.chains.back();

    gemmi::Residue res;
    res.name = "UNK";
    res.seqid = gemmi::SeqId(1, ' ');
    ch.residues.push_back(std::move(res));
    gemmi::Residue& r = ch.residues.back();

    for (const auto& a : m.atoms) {
        gemmi::Atom atom;
        atom.name = a.element;
        atom.element = gemmi::Element(a.element);
        atom.pos = gemmi::Position(a.atom_position.x(), a.atom_position.y(), a.atom_position.z());
        r.atoms.push_back(atom);
    }

    asc.mol = structure;
    asc.n_selected_atoms = static_cast<int>(m.atoms.size());
    asc.read_success = 1;

    return asc;
}

} // namespace coot