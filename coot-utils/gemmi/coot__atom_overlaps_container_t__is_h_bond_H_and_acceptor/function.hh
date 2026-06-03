#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include "coot/geometry/hb-types.hh"
#include "coot/geometry/protein-geometry.hh"

namespace coot {

class protein_geometry;

class atom_overlaps_container_t {
public:
  atom_overlaps_container_t(coot::protein_geometry* g = nullptr) : geom_p(g) {}

  coot::hb_t get_h_bond_type_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue);

  std::pair<bool, bool> is_h_bond_H_and_acceptor_gemmi(
      const gemmi::Atom& ligand_atom,
      const gemmi::Residue& ligand_res,
      const gemmi::Atom& env_atom,
      const gemmi::Residue& env_res
  );

private:
  coot::protein_geometry* geom_p;
};

inline coot::hb_t
atom_overlaps_container_t::get_h_bond_type_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue) {
   coot::hb_t type = coot::hb_t::HB_UNASSIGNED;
   if (geom_p) {
      std::string atom_name = atom.name;
      std::string res_name = residue.name;
      type = geom_p->get_h_bond_type(atom_name, res_name, coot::protein_geometry::IMOL_ENC_ANY);
   }
   return type;
}

} // namespace coot