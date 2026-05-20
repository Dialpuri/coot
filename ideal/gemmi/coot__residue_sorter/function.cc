#include "function.hh"
#include <string>

namespace coot {

bool residue_sorter_gemmi(const std::pair<bool, gemmi::CRA> &r1,
                          const std::pair<bool, gemmi::CRA> &r2) {
  const gemmi::Chain* chain1 = r1.second.chain;
  const gemmi::Chain* chain2 = r2.second.chain;
  const gemmi::Residue* res1 = r1.second.residue;
  const gemmi::Residue* res2 = r2.second.residue;

  if (!chain1 || !chain2 || !res1 || !res2) {
    return false;
  }

  std::string chain_id_1 = chain1->name;
  std::string chain_id_2 = chain2->name;
  if (chain_id_1 < chain_id_2) {
    return true;
  } else {
    if (chain_id_1 > chain_id_2) {
      return false;
    } else {
      int seq1 = res1->seqid.num.value;
      int seq2 = res2->seqid.num.value;
      if (seq1 < seq2) {
        return true;
      } else {
        if (seq1 > seq2) {
          return false;
        } else {
          char icode1 = res1->seqid.icode;
          char icode2 = res2->seqid.icode;
          if (icode1 < icode2) {
            return true;
          } else {
            if (icode1 > icode2) {
              return false;
            }
          }
        }
      }
    }
  }
  return false;
}

} // namespace coot