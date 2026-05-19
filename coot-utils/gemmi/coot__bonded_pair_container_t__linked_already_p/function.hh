#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

struct bonded_pair_t_gemmi {
    gemmi::Residue *res_1;
    gemmi::Residue *res_2;
    std::string link_type;
    bool is_fixed_first;
    bool is_fixed_second;
    
    bonded_pair_t_gemmi() : res_1(nullptr), res_2(nullptr), 
                            is_fixed_first(false), is_fixed_second(false) {}
    
    bonded_pair_t_gemmi(gemmi::Residue *r1, gemmi::Residue *r2, 
                        bool is_fixed_first_in, bool is_fixed_second_in,
                        const std::string &lt)
        : res_1(r1), res_2(r2), link_type(lt),
          is_fixed_first(is_fixed_first_in), is_fixed_second(is_fixed_second_in) {}
    
    bool matches(gemmi::Residue *r1, gemmi::Residue *r2) const {
        if (r1 == res_1 && r2 == res_2) {
            return true;
        } else if (r1 == res_2 && r2 == res_1) {
            return true;
        } else {
            return false;
        }
    }
};

class bonded_pair_container_t_gemmi {
public:
    std::vector<bonded_pair_t_gemmi> bonded_residues;
    
    unsigned int size() const { return bonded_residues.size(); }
    
    bool linked_already_p(gemmi::Residue *r1, gemmi::Residue *r2) const {
        for (size_t i = 0; i < bonded_residues.size(); i++) {
            if (((bonded_residues[i].res_1 == r1) &&
                 (bonded_residues[i].res_2 == r2)) ||
                ((bonded_residues[i].res_1 == r2) &&
                 (bonded_residues[i].res_2 == r1))) {
                return true;
            }
        }
        return false;
    }
};

} // namespace coot