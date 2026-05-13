#pragma once
#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Sort residues by sequence number (ascending), then by insertion code (descending)
inline std::vector<std::pair<gemmi::Residue*, int>>
sort_residues_by_seqno(gemmi::Residue** residues, int nResidues) {
    std::vector<std::pair<gemmi::Residue*, int>> v;
    
    // construct the vector
    for (int i = 0; i < nResidues; i++)
        v.push_back(std::make_pair(residues[i], i));
    
    // sort with comparison function (ascending by seqnum, descending by inscode)
    auto compare = [](const std::pair<gemmi::Residue*, int>& a,
                      const std::pair<gemmi::Residue*, int>& b) {
        int r1 = a.first->seqid.num.value;
        int r2 = b.first->seqid.num.value;
        
        if (r1 < r2) {
            return true;
        } else if (r1 > r2) {
            return false;
        } else {
            // Same seqnum: compare insertion code (descending order)
            char ins1 = a.first->seqid.icode;
            char ins2 = b.first->seqid.icode;
            return ins1 > ins2;
        }
    };
    
    std::sort(v.begin(), v.end(), compare);
    
    return v;
}

} // namespace util
} // namespace coot