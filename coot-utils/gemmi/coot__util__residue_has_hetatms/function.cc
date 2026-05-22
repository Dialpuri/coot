#include "function.hh"

namespace coot { namespace util {

int residue_has_hetatms_gemmi(const gemmi::Residue &res) {
    if (res.atoms.empty())
        return -1;
    if (res.het_flag == 'H')
        return 1;
    return 0;
}

}} // namespace coot::util