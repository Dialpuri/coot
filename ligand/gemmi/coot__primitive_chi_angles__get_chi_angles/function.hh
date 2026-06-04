#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

#include "clipper/core/clipper_util.h"

// Co-dependencies
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__primitive_chi_angles__get_quads/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__primitive_chi_angles__atom_to_co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__primitive_chi_angles__get_quads_using_altconfs/gemmi/function.hh"

namespace coot {

struct alt_confed_chi_angles {
    std::string alt_conf;
    std::vector<std::pair<int, float>> chi_angles;

    alt_confed_chi_angles() {}

    alt_confed_chi_angles(const std::string& conf,
                          const std::vector<std::pair<int, float>>& angles)
        : alt_conf(conf), chi_angles(angles) {}
};

} // namespace coot

namespace coot {
namespace primitive_chi_angles {

// ── chi-angle atom table types (self-contained, no extern) ──
// We use coot::primitive_chi_angles::atom_name_quad from get_quads_using_altconfs.hh
// (fields: names[4])

struct residue_named_chi_angle_atom_name_quad_set_t {
    std::string residue_type;
    std::vector<coot::primitive_chi_angles::atom_name_quad> name_quad;
    residue_named_chi_angle_atom_name_quad_set_t() {}
    residue_named_chi_angle_atom_name_quad_set_t(const std::string& rt)
        : residue_type(rt) {}
    void add_quad(const std::string& a1, const std::string& a2,
                  const std::string& a3, const std::string& a4) {
        name_quad.push_back(coot::primitive_chi_angles::atom_name_quad(a1, a2, a3, a4));
    }
};

// Global table — filled by setup_chi_atom_quads()
inline std::vector<residue_named_chi_angle_atom_name_quad_set_t>& chi_table() {
    static std::vector<residue_named_chi_angle_atom_name_quad_set_t> table;
    return table;
}

// Add a chi-atom quad to the table
inline void add_chi_quad_table(const std::string& residue_type,
                               const std::string& atom_name_1,
                               const std::string& atom_name_2,
                               const std::string& atom_name_3,
                               const std::string& atom_name_4) {
    auto& tbl = chi_table();
    for (auto& entry : tbl) {
        if (entry.residue_type == residue_type) {
            entry.add_quad(atom_name_1, atom_name_2, atom_name_3, atom_name_4);
            return;
        }
    }
    auto entry = residue_named_chi_angle_atom_name_quad_set_t(residue_type);
    entry.add_quad(atom_name_1, atom_name_2, atom_name_3, atom_name_4);
    tbl.push_back(entry);
}

// ── get_chi_angles_gemmi ────────────────────────────────────
std::vector<coot::alt_confed_chi_angles> get_chi_angles_gemmi(const gemmi::Residue& residue);

} // namespace primitive_chi_angles
} // namespace coot