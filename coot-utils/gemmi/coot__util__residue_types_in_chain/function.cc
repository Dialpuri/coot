#include "function.hh"
#include <set>
#include <algorithm>

namespace coot {
namespace util {

std::vector<std::string>
residue_types_in_chain_gemmi(const gemmi::Chain &chain) {
   std::set<std::string> unique_types;
   for (const auto &residue : chain.residues) {
      unique_types.insert(residue.name);
   }
   return std::vector<std::string>(unique_types.begin(), unique_types.end());
}

} // namespace util
} // namespace coot