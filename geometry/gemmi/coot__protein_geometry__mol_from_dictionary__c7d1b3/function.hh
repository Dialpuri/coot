#pragma once

#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__get_residue/gemmi/function.hh"

namespace coot {

// Gemmi port of coot::protein_geometry::mol_from_dictionary
// Builds a gemmi::Structure containing chain "A" with a single residue
// derived from the dictionary entry for comp_id.
// Returns a heap-allocated Structure (caller must delete) or nullptr.
gemmi::Structure* mol_from_dictionary_gemmi(
    protein_geometry& geom,
    const std::string& comp_id,
    int imol_enc,
    bool idealised_flag)
{
    gemmi::Structure* mol = nullptr;
    gemmi::Residue* residue_p = nullptr;
    float b_factor = 30.0f;

    residue_p = coot::get_residue_gemmi(
        geom, comp_id, imol_enc, idealised_flag, true, b_factor);

    if (residue_p) {
        gemmi::Chain chain;
        chain.name = "A";
        chain.residues.push_back(*residue_p);
        delete residue_p;

        gemmi::Model model;
        model.chains.push_back(chain);

        mol = new gemmi::Structure();
        mol->models.push_back(model);
    } else {
        std::cout << "WARNING:: Null residue in mol_from_dictionary() for comp_id "
                  << comp_id << std::endl;
    }

    std::cout << "DEBUG:: mol_from_dictionary() returns " << mol << std::endl;
    return mol;
}

} // namespace coot