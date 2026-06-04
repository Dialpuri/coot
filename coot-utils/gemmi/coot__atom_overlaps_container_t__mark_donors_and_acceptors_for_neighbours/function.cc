#include "function.hh"
#include <iostream>

namespace coot {

void mark_donors_and_acceptors_for_neighbours_gemmi(
    const protein_geometry* geom_p,
    const std::vector<gemmi::Residue*>& neighbours,
    std::unordered_map<const gemmi::Atom*, hb_t>& atom_hb_types,
    bool& molecule_has_hydrogens
) {

   for (unsigned int i = 0; i < neighbours.size(); i++) {
      try {
         std::string resname(neighbours[i]->name);
         std::pair<bool, dictionary_residue_restraints_t> d =
             geom_p->get_monomer_restraints(resname, protein_geometry::IMOL_ENC_ANY);
         if (!d.first) continue;

         const dictionary_residue_restraints_t& dict = d.second;

         for (gemmi::Atom& atom : neighbours[i]->atoms) {
            std::string atom_name(atom.name);
            std::string ele = atom.element.name();

            if (ele == "H") {
               molecule_has_hydrogens = true;
               std::string heavy_neighb_of_H_atom = dict.get_bonded_atom(atom_name);
               if (!heavy_neighb_of_H_atom.empty()) {
                  std::string neigh_energy_type = dict.type_energy(heavy_neighb_of_H_atom);
                  energy_lib_atom neighb_ela = geom_p->get_energy_lib_atom(neigh_energy_type);
                  hb_t neighb_hb_type = neighb_ela.hb_type;

                  if (neighb_hb_type == HB_DONOR) {
                     atom_hb_types[&atom] = HB_HYDROGEN;
                  }
                  if (neighb_hb_type == HB_BOTH) {
                     atom_hb_types[&atom] = HB_HYDROGEN;
                  }
                  if (neigh_energy_type == "NR5") {
                     atom_hb_types[&atom] = HB_HYDROGEN;
                  }
               } else {
                  std::cout << "ERROR:: mark_donors_and_acceptors_for_neighbours_gemmi() oops empty heavy for H-atom "
                            << atom_name << std::endl;
               }
            } else {
               std::string energy_type = dict.type_energy(atom_name);
               energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
               hb_t hb_type = ela.hb_type;

               atom_hb_types[&atom] = hb_type;

               if (energy_type == "NR5") {
                  atom_hb_types[&atom] = HB_DONOR;
               }
            }
         }
      }
      catch (const std::out_of_range& ex) {
         std::cout << "Opps " << ex.what() << std::endl;
      }
   }
}

} // namespace coot