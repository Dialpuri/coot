#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// Container for matching atom pairs (from first selection to second selection)
struct match_container_for_residues_t {
    struct AtomPair {
        gemmi::CRA first;
        gemmi::CRA second;
    };
    
    std::vector<AtomPair> atom_pairs;
    
    // Add a matching pair
    void add(const gemmi::CRA& first, const gemmi::CRA& second) {
        atom_pairs.push_back({first, second});
    }
    
    // Delete residues upstream (before) the matching atom in the selection
    // atoms0 = atoms from first selection
    // atoms1 = atoms from second selection
    // from_first: if true, use atoms0 and match against atom_pairs[ip].first
    //             if false, use atoms1 and match against atom_pairs[ip].second
    void delete_upstream_gemmi(gemmi::Structure& st, bool from_first,
                                const std::vector<gemmi::CRA>& atoms0,
                                const std::vector<gemmi::CRA>& atoms1);
};

} // namespace coot