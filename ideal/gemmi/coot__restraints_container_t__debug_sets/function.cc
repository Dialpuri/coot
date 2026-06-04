#include "function.hh"
#include "geometry/residue-and-atom-specs.hh"
#include <iostream>
#include <sstream>

namespace coot {

void restraints_container_t::debug_sets_gemmi() const {

   std::cout << "-------------------- in debug_sets() residues_vec: " << std::endl;
   for (std::size_t i=0; i<residues_vec.size(); i++)
      std::cout << "   " << residues_vec[i].first << " "
                << residue_spec_t(1, residues_vec[i].second ? residues_vec[i].second->name : "",
                                  residues_vec[i].second ? residues_vec[i].second->seqid.num.value : 0,
                                  residues_vec[i].second ? std::string(1, residues_vec[i].second->seqid.icode) : "")
                << std::endl;

   std::map<gemmi::Residue*, std::set<gemmi::Residue*>>::const_iterator it;
   for (it=fixed_neighbours_set.begin(); it!=fixed_neighbours_set.end(); it++) {
      gemmi::Residue* res_ptr = it->first;
      std::cout << "   Moving residue "
                << residue_spec_t(1,
                                  res_ptr ? res_ptr->name : "",
                                  res_ptr ? res_ptr->seqid.num.value : 0,
                                  res_ptr ? std::string(1, res_ptr->seqid.icode) : "")
                << std::endl;
      const std::set<gemmi::Residue*> &s = it->second;
      std::set<gemmi::Residue*>::const_iterator its;
      for (its=s.begin(); its!=s.end(); its++) {
         gemmi::Residue* n_ptr = *its;
         std::cout << "     fixed neigb: "
                   << residue_spec_t(1,
                                     n_ptr ? n_ptr->name : "",
                                     n_ptr ? n_ptr->seqid.num.value : 0,
                                     n_ptr ? std::string(1, n_ptr->seqid.icode) : "")
                   << std::endl;
      }
   }

   if (mol) {
      for (gemmi::Chain& chain_p : mol->chains) {
         std::cout << "   Chain " << chain_p.name << std::endl;
         for (gemmi::Residue& residue_p : chain_p.residues) {
            std::cout << "      "
                      << residue_spec_t(1,
                                        residue_p.name,
                                        residue_p.seqid.num.value,
                                        std::string(1, residue_p.seqid.icode))
                      << std::endl;
            for (gemmi::Atom& at : residue_p.atoms) {
               std::cout << "   "
                         << atom_spec_t(chain_p.name,
                                        residue_p.seqid.num.value,
                                        std::string(1, residue_p.seqid.icode),
                                        at.name,
                                        std::string(1, at.altloc))
                         << std::endl;
            }
         }
      }
   }

}

} // namespace coot