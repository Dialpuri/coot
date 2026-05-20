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
    
    // Delete residues downstream from matching atoms
    // atom_vectors[0] = atoms from first selection
    // atom_vectors[1] = atoms from second selection
    // from_first: if true, use first selection as matching reference
    //             if false, use second selection as matching reference
    void delete_downstream_gemmi(gemmi::Structure& st, bool from_first,
                                  const std::vector<gemmi::CRA>& atom_vectors_0,
                                  const std::vector<gemmi::CRA>& atom_vectors_1);
};

} // namespace coot