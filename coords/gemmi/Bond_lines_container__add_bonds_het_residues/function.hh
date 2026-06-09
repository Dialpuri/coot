#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include <set>
#include "coords/Bond_lines.hh"
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/mmdb-extras.hh"
#include "coot-utils/atom-selection-container.hh"
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__add_half_bonds/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__atom_colour/gemmi/function.hh"

// Check if a residue name is a standard amino acid or nucleotide
inline bool is_std_aa_or_nuc(const std::string &name) {
    static const char *std_aa[] = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
    };
    for (const char *n : std_aa) if (name == n) return true;
    static const char *std_nuc[] = {"A","G","C","U","T","DA","DG","DC","DU","DT"};
    for (const char *n : std_nuc) if (name == n) return true;
    return false;
}

inline bool is_hydrogen_element(const std::string &el) {
    gemmi::Element elem(el);
    return elem.is_hydrogen();
}

// gemmi port of Bond_lines_container::add_bonds_het_residues
// Takes Bond_lines_container_gemmi (for bonds storage + add_bond_half_gemmi)
// and explicit parameters for what the original accessed via private members.
void add_bonds_het_residues_gemmi(
    Bond_lines_container_gemmi &blc,
    const std::vector<std::pair<bool, gemmi::Residue *>> &het_residues,
    const atom_selection_container_t &atom_sel,
    int imol,
    int atom_colour_type,
    short int have_udd_handle,
    int udd_bond_handle,
    int udd_atom_index_handle,
    int udd_user_defined_atom_colour_index_handle,
    const coot::protein_geometry *geom,
    bool have_dictionary,
    bool do_bonds_to_hydrogens)
{
    if (!het_residues.empty()) {
        for (unsigned int ires = 0; ires < het_residues.size(); ires++) {
            if (het_residues[ires].first) {
                gemmi::Residue *res = het_residues[ires].second;
                std::string res_name = res->name;

                // Skip standard amino acids and nucleotides (original: isAminoacid/isNucleotide)
                if (is_std_aa_or_nuc(res_name)) continue;

                std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
                    geom->get_monomer_restraints_at_least_minimal(res_name, imol);

                if (!restraints.first) {
                    std::cout << "Oooppps!  No bonding rules for residue type :" << res_name
                              << ": missing bonds! " << std::endl;
                } else {
                    std::vector<gemmi::Atom *> residue_atoms;
                    for (auto &atom : res->atoms) {
                        residue_atoms.push_back(&atom);
                    }
                    int n_atoms = static_cast<int>(residue_atoms.size());
                    int model_number = 1;

                    for (unsigned int ib = 0; ib < restraints.second.bond_restraint.size(); ib++) {
                        std::string atom_name_1 = restraints.second.bond_restraint[ib].atom_id_1_4c();
                        std::string atom_name_2 = restraints.second.bond_restraint[ib].atom_id_2_4c();
                        std::string bt = restraints.second.bond_restraint[ib].type();

                        for (int iat = 0; iat < n_atoms; iat++) {
                            std::string residue_atom_name_1(residue_atoms[iat]->name);
                            if (atom_name_1 == residue_atom_name_1) {
                                for (int jat = 0; jat < n_atoms; jat++) {
                                    std::string residue_atom_name_2(residue_atoms[jat]->name);
                                    if (atom_name_2 == residue_atom_name_2) {
                                        char aloc_1 = residue_atoms[iat]->altloc;
                                        char aloc_2 = residue_atoms[jat]->altloc;
                                        if (aloc_1 == aloc_2 || aloc_1 == '\0' || aloc_2 == '\0') {
                                            coot::Cartesian p1(residue_atoms[iat]->pos.x,
                                                               residue_atoms[iat]->pos.y,
                                                               residue_atoms[iat]->pos.z);
                                            coot::Cartesian p2(residue_atoms[jat]->pos.x,
                                                               residue_atoms[jat]->pos.y,
                                                               residue_atoms[jat]->pos.z);

                                            int iat_1_atom_index = -1;
                                            int iat_2_atom_index = -1;

                                            std::string element_1 = residue_atoms[iat]->element.name();
                                            std::string element_2 = residue_atoms[jat]->element.name();

                                            gemmi::CRA cra1{nullptr, res, residue_atoms[iat]};
                                            gemmi::CRA cra2{nullptr, res, residue_atoms[jat]};

                                            if (element_1 != element_2) {
                                                if (!is_hydrogen_element(element_1) && !is_hydrogen_element(element_2)) {
                                                    if (bt == "double" || bt == "deloc" || bt == "triple") {
                                                        // add_double_bond / add_triple_bond — no compatible gemmi port
                                                    } else {
                                                        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                                        add_half_bonds_gemmi(
                                                            blc, p1, p2,
                                                            cra1, cra2,
                                                            cc, model_number,
                                                            iat_1_atom_index, iat_2_atom_index,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            nullptr, false, false);
                                                    }
                                                } else {
                                                    if (do_bonds_to_hydrogens) {
                                                        if (res_name == "HOH" || res_name == "DOD") {
                                                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                                            add_half_bonds_gemmi(
                                                                blc, p1, p2,
                                                                cra1, cra2,
                                                                cc, model_number,
                                                                iat_1_atom_index, iat_2_atom_index,
                                                                atom_colour_type,
                                                                udd_user_defined_atom_colour_index_handle,
                                                                nullptr, true, true);
                                                        } else {
                                                            blc.add_bond_half_gemmi(
                                                                HYDROGEN_GREY_BOND, p1, p2,
                                                                graphics_line_t::SINGLE,
                                                                model_number,
                                                                iat_1_atom_index,
                                                                iat_2_atom_index,
                                                                false, false);
                                                        }
                                                    }
                                                }
                                            } else {
                                                int col = atom_colour_gemmi(cra1,
                                                                            atom_colour_type,
                                                                            udd_user_defined_atom_colour_index_handle,
                                                                            nullptr);
                                                if (bt == "double" || bt == "deloc" || bt == "triple") {
                                                    // skip
                                                } else {
                                                    blc.add_bond_half_gemmi(
                                                        col, p1, p2,
                                                        graphics_line_t::SINGLE,
                                                        model_number,
                                                        iat_1_atom_index,
                                                        iat_2_atom_index,
                                                        false, false);
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
    }
}