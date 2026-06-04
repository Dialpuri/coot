#pragma once

#include <gemmi/model.hpp>
#include <gemmi/to_pdb.hpp>
#include <fstream>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__correct_link_distances/gemmi/function.hh"

namespace coot {
namespace util {

// No gemmi port — cis-peptide removal is an MMDB-specific concept
// (mmdb::CisPep).  Gemmi has no equivalent, so this is a no-op.
inline void remove_wrong_cis_peptides_gemmi(gemmi::Structure&) {
    // no-op: gemmi has no cis-peptide metadata to clean
}

} // namespace util
} // namespace coot

namespace coot {

inline int write_coords_pdb_gemmi(gemmi::Structure& st, const std::string& file_name) {
    coot::util::remove_wrong_cis_peptides_gemmi(st);
    coot::util::correct_link_distances_gemmi(st);

    std::ofstream ofs(file_name);
    if (!ofs) {
        return 1; // failure
    }
    gemmi::write_pdb(st, ofs);
    ofs.close();
    return 0;
}

} // namespace coot