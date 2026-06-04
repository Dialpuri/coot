#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace pli {

class pi_stacking_container_t {
public:
    // Gemmi version: takes a gemmi::Residue reference instead of mmdb::Residue*
    // Returns positions of cation atoms (NZ for LYS, NH1/NH2 for ARG)
    std::vector<clipper::Coord_orth>
    get_cation_atom_positions_gemmi(const gemmi::Residue &res) const {
        std::vector<clipper::Coord_orth> v;

        // Trim whitespace from atom name for format-independent comparison
        auto trim = [](const std::string& s) {
            size_t start = s.find_first_not_of(' ');
            if (start == std::string::npos) return std::string();
            size_t end = s.find_last_not_of(' ');
            return s.substr(start, end - start + 1);
        };

        std::string res_name = res.name;

        if (res_name == "LYS") {
            for (const gemmi::Atom &atom : res.atoms) {
                std::string atom_name = trim(atom.name);
                if (atom_name == "NZ") {
                    clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                    v.push_back(pt);
                }
            }
        }

        if (res_name == "ARG") {
            for (const gemmi::Atom &atom : res.atoms) {
                std::string atom_name = trim(atom.name);
                if ((atom_name == "NH1") || (atom_name == "NH2")) {
                    clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                    v.push_back(pt);
                }
            }
        }

        return v;
    }
};

} // namespace pli