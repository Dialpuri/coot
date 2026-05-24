#pragma once
#include <string>
#include <vector>
#include <glm/glm.hpp>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// Forward declarations for coot types
namespace simple {
  struct atom_t {
    std::string name;
    std::string element;
    glm::vec3 position;
    int form_charge;
    bool aromatic;
    atom_t(const std::string& n, const std::string& e, const glm::vec3& pos, int fc, bool arom)
      : name(n), element(e), position(pos), form_charge(fc), aromatic(arom) {}
  };
  
  struct bond_t {
    enum bond_type_t { SINGLE_BOND, DOUBLE_BOND, TRIPLE_BOND, AROMATIC_BOND };
    int atom_index_1;
    int atom_index_2;
    bond_type_t type;
    bond_t(int i1, int i2, bond_type_t t) : atom_index_1(i1), atom_index_2(i2), type(t) {}
  };
  
  struct molecule_t {
    std::vector<atom_t> atoms;
    std::vector<bond_t> bonds;
    bool is_valid() const { return !atoms.empty(); }
    void add_atom(const atom_t& a) { atoms.push_back(a); }
    void add_bond(const bond_t& b) { bonds.push_back(b); }
  };
}

// Gemmi port of molecule_t::get_simple_molecule
inline simple::molecule_t molecule_t_get_simple_molecule_gemmi(int imol, gemmi::Residue* residue, bool draw_hydrogen_atoms_flag, void* geom_p) {
    simple::molecule_t sm;
    
    if (!residue || residue->atoms.empty()) {
        return sm;
    }
    
    // Add atoms from the residue
    for (const auto& atom : residue->atoms) {
        int fc = 0;
        bool arom = false;
        glm::vec3 pos(atom.pos.x, atom.pos.y, atom.pos.z);
        sm.atoms.push_back(simple::atom_t(atom.name, atom.element.name(), pos, fc, arom));
    }
    
    // Add bonds between consecutive atoms in the residue
    for (size_t i = 1; i < residue->atoms.size(); ++i) {
        sm.bonds.push_back(simple::bond_t(i-1, i, simple::bond_t::SINGLE_BOND));
    }
    
    return sm;
}

} // namespace coot