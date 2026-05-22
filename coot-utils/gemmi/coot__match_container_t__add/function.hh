#pragma once
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class match_container_for_residues_t_gemmi {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_pairs;

    match_container_for_residues_t_gemmi()
        : residue_1(nullptr), residue_2(nullptr) {}

    match_container_for_residues_t_gemmi(gemmi::Residue *r1, gemmi::Residue *r2)
        : residue_1(r1), residue_2(r2) {}

    void add(gemmi::Atom *at_1, gemmi::Atom *at_2) {
        atom_pairs.emplace_back(at_1, at_2);
    }
};

class match_container_t {
public:
    std::vector<match_container_for_residues_t_gemmi> matches;

    void add_gemmi(gemmi::Atom *at_1, gemmi::Atom *at_2, const gemmi::Model &model) {
        gemmi::Residue *res_1 = nullptr;
        gemmi::Residue *res_2 = nullptr;

        // Search for parent residues in the model
        for (const gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : const_cast<gemmi::Chain&>(chain).residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (&atom == at_1) res_1 = &res;
                    if (&atom == at_2) res_2 = &res;
                }
            }
        }

        if (res_1 && res_2) {
            bool added = false;
            for (size_t i = 0; i < matches.size(); i++) {
                if (matches[i].residue_1 == res_1) {
                    if (matches[i].residue_2 == res_2) {
                        matches[i].add(at_1, at_2);
                        added = true;
                        break;
                    }
                }
            }
            if (!added) {
                // make a new one
                match_container_for_residues_t_gemmi m(res_1, res_2);
                m.add(at_1, at_2);
                matches.push_back(m);
            }
        }
    }
};

}  // namespace coot