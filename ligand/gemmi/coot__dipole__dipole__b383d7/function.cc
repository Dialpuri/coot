#include "function.hh"

namespace coot {
namespace dipole {

dipole_state_t dipole_gemmi(
    const coot::dictionary_residue_restraints_t &rest,
    gemmi::CRA cra) {

  std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> dict_res_pairs;
  dict_res_pairs.push_back(std::make_pair(rest, cra));

  dipole_state_t result;
  coot::dipole::init_gemmi(result, dict_res_pairs);

  return result;
}

} // namespace dipole
} // namespace coot
