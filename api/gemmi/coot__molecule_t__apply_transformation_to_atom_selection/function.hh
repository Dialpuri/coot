#pragma once
#include <string>
#include <vector>
#include <tuple>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>

namespace coot {
namespace molecule_t {

// Parse selection string like "//A" or "//A/10" and return atoms matching it
// Returns vector of (chain_name, residue_seq, atom_name) for matching atoms
std::vector<std::tuple<std::string, int, std::string>>
parse_selection(const std::string& selection_cid, const gemmi::Structure& st);

// Apply transformation to selected atoms
// selection_cid: selection string (e.g. "//A", "//A/10")
// n_atoms: expected number of atoms (for validation)
// rotation_centre: centre of transformation (gemmi::Vec3)
// rtop: transformation (gemmi::Transform)
// Returns: number of atoms moved
int apply_transformation_to_atom_selection_gemmi(
    const std::string& selection_cid,
    int n_atoms,
    const gemmi::Vec3& rotation_centre,
    const gemmi::Transform& rtop,
    gemmi::Structure& st);

} // namespace molecule_t
} // namespace coot