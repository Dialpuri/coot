#include "function.hh"
#include <gemmi/pdb.hpp>
#include <string>
#include <vector>

namespace coot {

// Helper to normalize insertion code: MMDB "" → gemmi ' '
static std::string norm_icode(const std::string& s) {
  return s.empty() ? std::string(" ") : s;
}

// Helper to trim spaces from atom name (gemmi stores trimmed names)
static std::string trim_atom_name(const std::string& name) {
  // Gemmi stores atom names without padding, so we just need to match as-is
  // The spec might have padding like " CA ", so trim it for comparison
  size_t start = name.find_first_not_of(" ");
  size_t end = name.find_last_not_of(" ");
  if (start == std::string::npos) return "";
  return name.substr(start, end - start + 1);
}

// matches_spec gemmi port - uses atom_spec_t::matches_spec_gemmi
static bool matches_spec_gemmi(const atom_spec_t& spec, const gemmi::CRA& cra) {
  if (!cra.residue || !cra.atom)
    return false;
  
  // Chain ID match (chain_id in spec maps to parent chain name)
  if (!spec.chain_id.empty() && cra.chain && cra.chain->name != spec.chain_id)
    return false;
  
  // Residue number match
  if (cra.residue->seqid.num.value != spec.res_no)
    return false;
  
  // Insertion code match
  std::string res_icode_str(1, cra.residue->seqid.icode);
  if (norm_icode(spec.ins_code) != res_icode_str)
    return false;
  
  // Atom name match - gemmi stores atom names without padding
  // Trim spaces from spec.atom_name for comparison
  std::string spec_name = trim_atom_name(spec.atom_name);
  if (cra.atom->name != spec_name)
    return false;
  
  // Alt loc match
  if (!spec.alt_conf.empty() && cra.atom->altloc != spec.alt_conf[0])
    return false;
  
  return true;
}

namespace torsion_general {

int atom_index_gemmi(const atom_spec_t& spec, gemmi::Residue& residue) {
  const int FAIL = -1;
  
  for (size_t i = 0; i < residue.atoms.size(); ++i) {
    // Build CRA with nullptr for chain (we don't have parent context here)
    gemmi::CRA cra{nullptr, &residue, &const_cast<gemmi::Atom&>(residue.atoms[i])};
    if (matches_spec_gemmi(spec, cra)) {
      return static_cast<int>(i);
    }
  }
  return FAIL;
}

} // namespace torsion_general
} // namespace coot