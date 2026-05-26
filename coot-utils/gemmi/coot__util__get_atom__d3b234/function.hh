#pragma once
#include <string>
#include <gemmi/model.hpp>
#include "coot/geometry/residue-and-atom-specs.hh"

namespace coot {
namespace util {

// Ported from MMDB to gemmi: searches for atom within a given residue
// Returns nullptr if not found or if residue is null
inline gemmi::Atom *get_atom_gemmi(const atom_spec_t &spec, gemmi::Residue *res) {
    if (!res) {
        return nullptr;
    }
    
    for (gemmi::Atom& test_at : res->atoms) {
        if (spec.atom_name == test_at.name) {
            std::string test_altconf(1, test_at.altloc);
            if (spec.alt_conf == test_altconf) {
                return &test_at;
            }
        }
    }
    return nullptr;
}

} // namespace util
} // namespace coot