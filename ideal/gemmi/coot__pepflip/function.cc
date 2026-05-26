#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__pepflip_standard/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__pepflip_internal_to_residue/gemmi/function.hh"

namespace coot {

int pepflip_gemmi(gemmi::Structure &st,
                  const std::string &chain_id,
                  int resno,
                  const std::string &ins_code,
                  const std::string &altconf) {
   int status = pepflip_standard_gemmi(st, chain_id, resno, ins_code, altconf);
   if (status)
      return status;
   else
      return pepflip_internal_to_residue_gemmi(st, chain_id, resno, ins_code, altconf);
}

} // namespace coot