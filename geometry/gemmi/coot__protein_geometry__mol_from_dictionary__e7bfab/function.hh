#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <memory>
#include <iostream>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__get_residue/gemmi/function.hh"

namespace coot {

std::unique_ptr<gemmi::Structure> mol_from_dictionary_gemmi(
    coot::protein_geometry& geom,
    const std::string &three_letter_code,
    int imol_enc,
    bool idealised_flag) {

    std::unique_ptr<gemmi::Structure> mol = nullptr;

    gemmi::Residue *residue_p = get_residue_gemmi(
        geom, three_letter_code, imol_enc, idealised_flag, true, 0.f);

    if (residue_p) {
        mol = std::make_unique<gemmi::Structure>();
        auto& model = mol->models.emplace_back();
        auto& chain = model.chains.emplace_back();
        chain.name = "A";
        chain.residues.push_back(*residue_p);
        delete residue_p;  // get_residue_gemmi returns a heap-allocated residue
    } else {
        std::cout << "WARNING:: protein_geometry::mol_from_dictionary(): "
                  << "Null residue in mol_from_dictionary() for "
                  << three_letter_code << std::endl;
    }

    return mol;
}

}  // namespace coot