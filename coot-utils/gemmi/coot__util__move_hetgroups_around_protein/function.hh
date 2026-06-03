#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <iostream>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__move_waters_around_protein/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_cell_symm/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__residue_has_hetatms/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_centre/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__symmetry_move_atoms/gemmi/function.hh"

namespace coot {
namespace util {

void move_hetgroups_around_protein_gemmi(gemmi::Structure& st) {

    // do the waters first
    coot::util::move_waters_around_protein_gemmi(st);

    try {
        // Now clipperize the variables.
        std::pair<clipper::Cell, clipper::Spacegroup> csp = coot::util::get_cell_symm_gemmi(st);
        clipper::Cell cell = csp.first;
        clipper::Spacegroup spacegroup = csp.second;

        if (cell.is_null()) {
            std::cout << "WARNING:: null cell in move_waters_around_protein" << std::endl;
        } else {
            if (spacegroup.is_null()) {
                std::cout << "WARNING:: null spgr in move_waters_around_protein" << std::endl;
            } else {

                for (gemmi::Model& model : st.models) {
                    // first find the protein coords
                    std::vector<clipper::Coord_orth> protein_coords;

                    for (gemmi::Chain& chain : model.chains) {
                        for (gemmi::Residue& residue : chain.residues) {
                            for (gemmi::Atom& atom : residue.atoms) {
                                if (residue.het_flag != 'H') {
                                    std::string element = atom.element.name();
                                    if (element != "C") {
                                        clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                                        protein_coords.push_back(pt);
                                    }
                                }
                            }
                        }
                    }
                    // OK, protein_coords is set.

                    // run over chains of the existing mol
                    for (gemmi::Chain& chain : model.chains) {
                        for (gemmi::Residue& residue : chain.residues) {
                            std::string residue_name(residue.name);
                            if (residue_name == "WAT" || residue_name == "HOH") {
                                // Waters are handled above.
                            } else {
                                if (coot::util::residue_has_hetatms_gemmi(residue)) {
                                    std::pair<bool, clipper::Coord_orth> centre =
                                        coot::util::get_residue_centre_gemmi(residue);
                                    if (centre.first) {

                                        // fill hetgroup_atoms:
                                        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>> hetgroup_atoms;
                                        for (gemmi::Atom& atom : residue.atoms) {
                                            clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
                                            hetgroup_atoms.emplace_back(&atom, co);
                                        }

                                        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>> atoms_moved =
                                            coot::util::symmetry_move_atoms_gemmi(
                                                protein_coords, hetgroup_atoms, cell, spacegroup);

                                        for (const auto& p : atoms_moved) {
                                            if (p.first) {
                                                p.first->pos.x = p.second.x();
                                                p.first->pos.y = p.second.y();
                                                p.first->pos.z = p.second.z();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch (const std::runtime_error &rte) {
        std::cout << rte.what() << std::endl;
    }
}

} // namespace util
} // namespace coot