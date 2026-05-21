#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Compare two residues by sequence number, then insertion code
// Returns true if a < b (a comes before b in sorting order)
inline bool compare_residues_gemmi(const std::pair<const gemmi::Residue*, int>& a,
                                   const std::pair<const gemmi::Residue*, int>& b) {
    int r1 = a.first->seqid.num.value;
    int r2 = b.first->seqid.num.value;

    if (r1 < r2) {
        return true;
    } else {
        if (r1 > r2) {
            return false;
        } else {
            char ins1 = a.first->seqid.icode;
            char ins2 = b.first->seqid.icode;
            if (ins1 > ins2) {
                return true;
            } else {
                return false;
            }
        }
    }
}

} // namespace util
} // namespace coot