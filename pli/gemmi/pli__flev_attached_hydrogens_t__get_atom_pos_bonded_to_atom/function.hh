#pragma once

#include "geometry/protein-geometry.hh"
#include "coot/coot-utils/coot-coord-extras.hh"
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <string>
#include <stdexcept>

namespace pli {

// Minimal gemmi-only declaration of flev_attached_hydrogens_t (avoids redefining the original in pli/flev-annotations.hh)
class flev_attached_hydrogens_t {
public:
    coot::dictionary_residue_restraints_t dict_res;

    flev_attached_hydrogens_t(coot::dictionary_residue_restraints_t& r) : dict_res(r) {}

    clipper::Coord_orth get_atom_pos_bonded_to_atom_gemmi(
        const gemmi::Atom *lig_at,
        const gemmi::Atom *H_at,
        const gemmi::Residue *ligand_residue,
        const coot::protein_geometry &geom) const;
};

} // namespace pli

clipper::Coord_orth
pli::flev_attached_hydrogens_t::get_atom_pos_bonded_to_atom_gemmi(
    const gemmi::Atom *lig_at,
    const gemmi::Atom *H_at,
    const gemmi::Residue *ligand_residue,
    const coot::protein_geometry &geom) const
{
    int imol = 0; // FIXME needs checking
    std::string res_name(ligand_residue->name);
    std::pair<bool, coot::dictionary_residue_restraints_t> p =
        geom.get_monomer_restraints_at_least_minimal(res_name, imol);

    if (!p.first) {
        std::string m = "No monomer type ";
        m += res_name;
        m += " found in dictionary";
        throw std::runtime_error(m);
    }

    const gemmi::Atom *bonded_atom = nullptr;
    std::string bonded_atom_name;
    // Pad gemmi atom names to 4-char format to match dictionary atom_id_1_4c() output
    std::string lig_at_name = coot::atom_id_mmdb_expand(lig_at->name);
    std::string H_at_name = coot::atom_id_mmdb_expand(H_at->name);

    for (unsigned int ibond = 0; ibond < p.second.bond_restraint.size(); ibond++) {
        std::string atom_name_1 = p.second.bond_restraint[ibond].atom_id_1_4c();
        std::string atom_name_2 = p.second.bond_restraint[ibond].atom_id_2_4c();
        if (atom_name_1 == lig_at_name) {
            if (atom_name_2 != H_at_name) {
                bonded_atom_name = atom_name_2;
                break;
            }
        }
        if (atom_name_2 == lig_at_name) {
            if (atom_name_1 != H_at_name) {
                bonded_atom_name = atom_name_1;
                break;
            }
        }
    }

    if (!bonded_atom_name.empty()) {
        for (const gemmi::Atom &atom : ligand_residue->atoms) {
            if (coot::atom_id_mmdb_expand(atom.name) == bonded_atom_name) {
                bonded_atom = &atom;
                break;
            }
        }
    }

    if (!bonded_atom) {
        std::string m = "No atom bonded to ";
        m += lig_at_name;
        m += " found in dictionary for ";
        m += res_name;
        throw std::runtime_error(m);
    }

    return clipper::Coord_orth(bonded_atom->pos.x, bonded_atom->pos.y, bonded_atom->pos.z);
}