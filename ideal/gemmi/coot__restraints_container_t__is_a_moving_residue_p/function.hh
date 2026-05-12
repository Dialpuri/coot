#pragma once
#include <set>
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

// Helper function to get a unique identifier for a gemmi::Residue
// to match the behavior of the original MMDB function that used pointers
struct residue_identifier {
    int seqnum;
    char icode;
    
    bool operator<(const residue_identifier& other) const {
        if (seqnum != other.seqnum) return seqnum < other.seqnum;
        return icode < other.icode;
    }
};

class restraints_container_t {
public:
    // The set of moving residues identified by their unique properties
    std::set<residue_identifier> residues_vec_moving_set;
    
    // Original function ported - checks if a residue is in the moving set
    bool is_a_moving_residue_p_gemmi(const gemmi::Residue* r) const {
        if (!r) return false;
        residue_identifier id;
        id.seqnum = r->seqid.num.value;
        id.icode = r->seqid.icode;
        return residues_vec_moving_set.find(id) != residues_vec_moving_set.end();
    }
};

} // namespace coot