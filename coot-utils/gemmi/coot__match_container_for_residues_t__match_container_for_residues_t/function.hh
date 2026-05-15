#pragma once
#include <utility>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class match_container_for_residues_t_gemmi {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::vector<gemmi::Residue *> fragment_1_res_vec;
    std::vector<gemmi::Residue *> fragment_2_res_vec;
    std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_pairs;

    match_container_for_residues_t_gemmi()
        : residue_1(nullptr), residue_2(nullptr) {}

    match_container_for_residues_t_gemmi(gemmi::Residue *r1, gemmi::Residue *r2)
        : residue_1(r1), residue_2(r2) {}

    void add(gemmi::Atom *at_1, gemmi::Atom *at_2) {
        atom_pairs.emplace_back(at_1, at_2);
    }
};

}  // namespace coot