#pragma once
#include <vector>
#include <stdexcept>
#include <string>
#include <gemmi/model.hpp>

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

        // Validate: all 4 residues must be distinct
        if (&res_0 == &res_1 || &res_0 == &res_2 || &res_0 == &res_3 ||
            &res_1 == &res_2 || &res_1 == &res_3 || &res_2 == &res_3) {
            throw(std::runtime_error("duplicate residues"));
        }

        // Find atoms in res_0 (should have C)
        for (auto& at : res_0.atoms) {
            if (at.name == "C") {
                v[0] = &at;
                break;
            }
        }

        // Find atoms in res_1 (should have N, C, O, CA)
        for (auto& at : res_1.atoms) {
            if (at.name == "N") v[1] = &at;
            if (at.name == "C") v[2] = &at;
            if (at.name == "O") v[3] = &at;
            if (at.name == "CA") ca_1 = &at;
        }

        // Find atoms in res_2 (should have N, H, C, CA)
        for (auto& at : res_2.atoms) {
            if (at.name == "N") v[4] = &at;
            if (at.name == "H") v[5] = &at;
            if (at.name == "C") v[6] = &at;
            if (at.name == "CA") ca_2 = &at;
        }

        // Find atoms in res_3 (should have N)
        for (auto& at : res_3.atoms) {
            if (at.name == "N") v[7] = &at;
        }

        if (!ca_1) throw(std::runtime_error("missing ca_1"));
        if (!ca_2) throw(std::runtime_error("missing ca_2"));

        // Check that we have all the atoms (except HN)
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
};

} // namespace coot