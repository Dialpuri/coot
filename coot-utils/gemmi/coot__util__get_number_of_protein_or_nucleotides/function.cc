#include "function.hh"
#include "utils/coot-utils.hh"

namespace coot {
namespace util {

std::pair<unsigned int, unsigned int>
get_number_of_protein_or_nucleotides_gemmi(gemmi::Chain* chain_p) {

   std::pair<unsigned int, unsigned int> n(0, 0);
   if (chain_p) {
      for (gemmi::Residue& residue : chain_p->residues) {
         std::string res_name(residue.name);
         if (is_standard_amino_acid_name(res_name))
            n.first++;
         if (is_standard_nucleotide_name(res_name))
            n.second++;
      }
   }

   return n;
}

} // namespace util
} // namespace coot