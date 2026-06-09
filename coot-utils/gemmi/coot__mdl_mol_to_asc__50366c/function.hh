#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "lidia-core/lbg-molfile.hh"

namespace coot {

gemmi::Structure mdl_mol_to_asc_gemmi(const lig_build::molfile_molecule_t& m) {
    gemmi::Structure st;
    if (m.atoms.empty()) {
        return st;
    }

    gemmi::Model model;
    model.num = 1;
    auto& chain = model.chains.emplace_back();
    chain.name = "A";

    gemmi::ResidueId rid;
    rid.seqid = gemmi::SeqId(1, ' ');
    rid.name = "UNL";
    gemmi::Residue res(rid);

    for (size_t i = 0; i < m.atoms.size(); i++) {
        const auto& a = m.atoms[i];
        gemmi::Atom atom;
        atom.element = gemmi::Element(a.element);
        atom.name = a.name;
        atom.pos = gemmi::Position(a.atom_position.x(), a.atom_position.y(), a.atom_position.z());
        atom.serial = static_cast<int>(i + 1);
        res.atoms.push_back(atom);
    }

    chain.residues.push_back(std::move(res));
    st.models.push_back(std::move(model));

    return st;
}

} // namespace coot
