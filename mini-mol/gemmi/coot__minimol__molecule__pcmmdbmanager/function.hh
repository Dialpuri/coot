#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace minimol {

struct atom {
    gemmi::Vec3 pos;
    double occupancy;
    double temperature_factor;
    std::string name;
    std::string element;
    char altloc;
    int int_user_data;
};

struct residue {
    std::string name;
    int seqnum;
    char ins_code;
    std::vector<atom> atoms;
};

struct fragment {
    std::string fragment_id;
    std::vector<residue> residues;
    
    int min_res_no() const {
        if (residues.empty()) return 0;
        int min_val = residues[0].seqnum;
        for (size_t i = 1; i < residues.size(); ++i) {
            if (residues[i].seqnum < min_val)
                min_val = residues[i].seqnum;
        }
        return min_val;
    }
    
    int max_residue_number() const {
        if (residues.empty()) return 0;
        int max_val = residues[0].seqnum;
        for (size_t i = 1; i < residues.size(); ++i) {
            if (residues[i].seqnum > max_val)
                max_val = residues[i].seqnum;
        }
        return max_val;
    }
    
    residue& operator[](int seqnum) {
        for (size_t i = 0; i < residues.size(); ++i) {
            if (residues[i].seqnum == seqnum)
                return residues[i];
        }
        return residues[0];
    }
    
    const residue& operator[](int seqnum) const {
        for (size_t i = 0; i < residues.size(); ++i) {
            if (residues[i].seqnum == seqnum)
                return residues[i];
        }
        return residues[0];
    }
};

struct molecule {
    std::vector<fragment> fragments;
    bool have_cell;
    bool have_spacegroup;
    double mmdb_cell[6];
    std::string mmdb_spacegroup;
    
    gemmi::Structure pcmmdbmanager_gemmi() const;
    
    fragment& operator[](unsigned int i) { return fragments[i]; }
    const fragment& operator[](unsigned int i) const { return fragments[i]; }
};

} // namespace minimol
} // namespace coot