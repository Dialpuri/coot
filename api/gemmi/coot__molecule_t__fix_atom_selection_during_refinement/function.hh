#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct atom_spec_t {
    std::string chain_name;
    int residue_number;
    std::string atom_name;
    char insertion_code;
    std::string element;
};

class molecule_t {
public:
    molecule_t(const gemmi::Structure& st) : structure_(st) {}
    
    void fix_atom_selection_during_refinement_gemmi(const std::string &atom_selection_cid);
    std::vector<atom_spec_t> get_fixed_atoms() const;
    
private:
    gemmi::Structure structure_;
};

} // namespace coot