#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "coot-utils/coot-coord-utils.hh"

namespace coot {

class crankshaft_set {
public:
    std::vector<gemmi::Atom*> v;
    gemmi::Atom* ca_1;
    gemmi::Atom* ca_2;

    crankshaft_set(gemmi::Residue& res_0, gemmi::Residue& res_1,
                   gemmi::Residue& res_2, gemmi::Residue& res_3) {
        v.resize(8, nullptr);
        ca_1 = nullptr;
        ca_2 = nullptr;

        if (&res_0 == &res_1 || &res_0 == &res_2 || &res_0 == &res_3 ||
            &res_1 == &res_2 || &res_1 == &res_3 || &res_2 == &res_3) {
            throw(std::runtime_error("duplicate residues"));
        }

        for (auto& at : res_0.atoms) {
            if (at.name == "C") {
                v[0] = &at;
                break;
            }
        }

        for (auto& at : res_1.atoms) {
            if (at.name == "N") v[1] = &at;
            if (at.name == "C") v[2] = &at;
            if (at.name == "O") v[3] = &at;
            if (at.name == "CA") ca_1 = &at;
        }

        for (auto& at : res_2.atoms) {
            if (at.name == "N") v[4] = &at;
            if (at.name == "H") v[5] = &at;
            if (at.name == "C") v[6] = &at;
            if (at.name == "CA") ca_2 = &at;
        }

        for (auto& at : res_3.atoms) {
            if (at.name == "N") v[7] = &at;
        }

        if (!ca_1) throw(std::runtime_error("missing ca_1"));
        if (!ca_2) throw(std::runtime_error("missing ca_2"));

        int n_atoms = 0;
        for (std::size_t i = 0; i < v.size(); i++) {
            if (v[i])
                n_atoms++;
        }
        if (n_atoms < 7) {
            std::string msg = "not enough atoms in crankshaft_set (";
            msg += std::to_string(n_atoms);
            msg += " < 7)";
            throw(std::runtime_error(msg));
        }
    }

    void move_the_atoms_gemmi(float ang) {
        if (v.size() > 5) {
            int indices[] = { 2, 3, 4, 5 };
            clipper::Coord_orth p_ca_1 = co_gemmi(ca_1);
            clipper::Coord_orth p_ca_2 = co_gemmi(ca_2);
            clipper::Coord_orth dir = p_ca_2 - p_ca_1;
            for (std::size_t i = 0; i < 4; i++) {
                gemmi::Atom* at = v[indices[i]];
                if (at) {
                    clipper::Coord_orth at_pos = co_gemmi(at);
                    clipper::Coord_orth at_pos_new = coot::util::rotate_around_vector(dir, at_pos, p_ca_1, ang);
                    at->pos.x = at_pos_new.x();
                    at->pos.y = at_pos_new.y();
                    at->pos.z = at_pos_new.z();
                }
            }
        }
    }
};

} // namespace coot