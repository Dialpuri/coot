#include "function.hh"
#include <iostream>

namespace coot {

short int progressive_residues_in_chain_check_gemmi(const gemmi::Chain& chain) {
  int previous_seq_num = -9999;

  for (const gemmi::Residue& res : chain.residues) {
    int this_seq_no = res.seqid.num.value;
    if (!(this_seq_no >= (previous_seq_num + 1))) {
      return 0;
    }
    previous_seq_num = this_seq_no;
  }

  return 1;
}

}  // namespace coot