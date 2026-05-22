#include "function.hh"

namespace coot {
namespace util {

gemmi::Residue *
previous_residue_gemmi(const gemmi::CRA& cra) {
  if (!cra.residue || !cra.chain)
    return nullptr;

  const std::vector<gemmi::Residue>& residues = cra.chain->residues;
  auto it = std::find(residues.begin(), residues.end(), *cra.residue);
  if (it == residues.end() || it == residues.begin())
    return nullptr;

  return const_cast<gemmi::Residue*>(&*(it - 1));
}

} // namespace util
} // namespace coot