#pragma once

#include <gemmi/model.hpp>
#include "geometry/hb-types.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

class protein_geometry;

class atom_overlaps_container_t {
  public:
    atom_overlaps_container_t(coot::protein_geometry* g) : geom_p(g) {}
    coot::hb_t get_h_bond_type_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue);

  private:
    coot::protein_geometry* geom_p;
};

inline coot::hb_t
atom_overlaps_container_t::get_h_bond_type_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue) {

   coot::hb_t type = coot::hb_t::HB_UNASSIGNED;
   std::string atom_name = atom.name;
   std::string res_name = residue.name;
   type = geom_p->get_h_bond_type(atom_name, res_name, coot::protein_geometry::IMOL_ENC_ANY); // heavyweight

   return type;
}

} // namespace coot