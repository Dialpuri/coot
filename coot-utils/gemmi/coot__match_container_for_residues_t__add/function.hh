#pragma once
#include <vector>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

class match_container_for_residues_t {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_pairs;

    match_container_for_residues_t(gemmi::Residue *r1, gemmi::Residue *r2)
        : residue_1(r1), residue_2(r2) {}
    match_container_for_residues_t() : residue_1(nullptr), residue_2(nullptr) {}

    void add_gemmi(gemmi::Atom *at_1, gemmi::Atom *at_2) {
        std::pair<gemmi::Atom *, gemmi::Atom *> p(at_1, at_2);
        atom_pairs.push_back(p);
    }
};

}  // namespace coot