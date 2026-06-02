#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__residue_names_with_no_dictionary/gemmi/function.hh"

namespace coot {

bool read_extra_dictionaries_for_molecule_gemmi(
    protein_geometry& geom,
    const gemmi::Structure* st,
    int imol_no,
    int* read_number_p) {

    if (!st) return false;

    std::vector<std::string> v = residue_names_with_no_dictionary_gemmi(geom, st, imol_no);

    bool success = true;
    for (std::size_t i = 0; i < v.size(); i++) {
        const std::string& rn = v[i];
        int success_for_residue = geom.try_dynamic_add(rn, *read_number_p);
        if (success_for_residue == 0)
            success = false;
        *read_number_p += 1;
    }

    return success;
}

} // namespace coot