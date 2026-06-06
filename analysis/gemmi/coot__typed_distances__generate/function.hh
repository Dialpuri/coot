#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__typed_distances__find_residues_within_window/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__typed_distances__get_type/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__distance/gemmi/function.hh"

namespace coot {
namespace typed_distances {

// get_atom_pair_bin_id inline helper (no port yet)
inline int get_atom_pair_bin_id(const atom_type_t& t1, const atom_type_t& t2) {
    int bin_id = -1;
    if (t1 == atom_type_t(C)) {
        if (t2 == atom_type_t(C)) bin_id = 0;
        if (t2 == atom_type_t(O)) bin_id = 1;
        if (t2 == atom_type_t(N)) bin_id = 2;
    }
    if (t1 == atom_type_t(O)) {
        if (t2 == atom_type_t(C)) bin_id = 1;
        if (t2 == atom_type_t(O)) bin_id = 3;
        if (t2 == atom_type_t(N)) bin_id = 4;
    }
    if (t1 == atom_type_t(N)) {
        if (t2 == atom_type_t(C)) bin_id = 2;
        if (t2 == atom_type_t(O)) bin_id = 4;
        if (t2 == atom_type_t(N)) bin_id = 5;
    }
    return bin_id;
}

void generate_gemmi(
    gemmi::Model& model,
    const std::vector<gemmi::CRA>& selected_atoms,
    std::map<gemmi::CRA, std::map<int, std::vector<float>>, cra_less>& residue_distances_map,
    std::map<gemmi::CRA, std::vector<gemmi::CRA>, cra_less>& residues_within_window
);

}  // namespace typed_distances
}  // namespace coot