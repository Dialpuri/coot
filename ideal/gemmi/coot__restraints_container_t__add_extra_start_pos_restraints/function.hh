#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

// Forward declare coot types used in signatures
// These come from coot headers; we only need their declarations.

namespace coot {

struct atom_spec_t {
    int model_number;
    std::string chain_id;
    int res_no;
    std::string ins_code;
    std::string atom_name;
    std::string alt_conf;
};

struct extra_restraints_t {
    struct extra_start_pos_restraint_t {
        atom_spec_t atom_1;
        double esd;
        extra_start_pos_restraint_t(const atom_spec_t& a, double e)
            : atom_1(a), esd(e) {}
    };
    std::vector<extra_start_pos_restraint_t> start_pos_restraints;
};

// gemmi port of restraints_container_t::add_extra_start_pos_restraints
// Uses gemmi Model for atom/residue lookup instead of MMDB.
// Returns the number of start_pos restraints that would be added
// (matching non-fixed atoms found in the model).
//
// The gemmi version skips UDData/fixed_check (MMDB-only concepts) and
// instead returns matching information. Callers can use the returned
// vector to decide which restraints to add.

struct start_pos_match_result {
    gemmi::CRA cra;       // matched atom location (may be null pointers if partial match)
    int restraint_index;  // index into extra_restraints.start_pos_restraints
    double esd;           // from the original restraint spec
};

std::vector<start_pos_match_result>
add_extra_start_pos_restraints_gemmi(
    gemmi::Model& model,
    const extra_restraints_t& extra_restraints);

} // namespace coot