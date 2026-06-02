#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>

#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__atoms_match_dictionary__c725f6/gemmi/function.hh"

bool add_bond_by_dictionary_maybe_gemmi(
    int imol,
    const gemmi::Atom* atom_p_1,
    const gemmi::Atom* atom_p_2,
    gemmi::Residue* residue_1,
    gemmi::Residue* residue_2,
    std::vector<std::pair<bool, gemmi::Residue *>> *het_residues,
    coot::protein_geometry* geom,
    bool have_dictionary,
    bool atom1_is_het,
    bool atom2_is_het)
{
    bool bond_het_residue_by_dictionary = false;
    if (have_dictionary && geom)
        if (residue_1 == residue_2)
            if (atom1_is_het)
                if (atom2_is_het) {

                    gemmi::Residue* shared_residue = residue_1;

                    std::pair<bool, gemmi::Residue*> tp0(false, shared_residue);
                    std::pair<bool, gemmi::Residue*> tp1(true, shared_residue);

                    std::vector<std::pair<bool, gemmi::Residue *>>::const_iterator it_1 =
                        std::find(het_residues->begin(), het_residues->end(), tp0);

                    if (it_1 == het_residues->end()) {

                        std::vector<std::pair<bool, gemmi::Residue *>>::const_iterator it_2 =
                            std::find(het_residues->begin(), het_residues->end(), tp1);

                        if (it_2 == het_residues->end()) {

                            auto dict_result = geom->get_monomer_restraints(shared_residue->name, imol);

                            if (dict_result.first) {

                                if (coot::atoms_match_dictionary_gemmi(*shared_residue, true, true).first) {

                                    het_residues->push_back(tp1);
                                    bond_het_residue_by_dictionary = true;
                                } else {
                                    het_residues->push_back(tp0);
                                }
                            } else {
                                het_residues->push_back(tp0);
                            }
                        } else {
                            bond_het_residue_by_dictionary = true;
                        }
                    } else {
                        bond_het_residue_by_dictionary = false;
                    }
                }

    return bond_het_residue_by_dictionary;
}