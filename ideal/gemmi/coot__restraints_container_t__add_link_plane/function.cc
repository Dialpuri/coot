#include "function.hh"
#include <iostream>

namespace coot {

int add_link_plane_gemmi(
    std::string link_type,
    gemmi::CRA first,
    gemmi::CRA second,
    short int is_fixed_first_res,
    short int is_fixed_second_res,
    const coot::protein_geometry& geom,
    std::vector<gemmi_plane_restraint_data_t>& planes_out)
{
   bool debug = false;
   if (debug && first.chain && second.chain)
      std::cout << "DEBUG:: add_link_plane() ::::::::  for type " << link_type << " "
                << first.chain->name << " " << first.residue->seqid.num.value
                << " :" << first.residue->seqid.icode << ":"
                << " -> " << second.chain->name << " " << second.residue->seqid.num.value
                << " :" << second.residue->seqid.icode << ":" << std::endl;

   int n_plane = 0;

   int n_first_res_atoms = (first.residue ? (int)first.residue->atoms.size() : 0);
   int n_second_res_atoms = (second.residue ? (int)second.residue->atoms.size() : 0);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   // Build atom_indices_map: alt_loc -> vector of atom indices
   std::map<std::string, std::vector<int> > atom_indices_map;

   if (first.residue) {
      for (const auto& atom : first.residue->atoms) {
         std::string alt_loc;
         if (atom.altloc != '\0') {
            alt_loc = std::string(1, atom.altloc);
         }
         std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
         if (it == atom_indices_map.end()) {
            std::vector<int> v;
            atom_indices_map[alt_loc] = v;
         }
      }
   }

   if (second.residue) {
      for (const auto& atom : second.residue->atoms) {
         std::string alt_loc;
         if (atom.altloc != '\0') {
            alt_loc = std::string(1, atom.altloc);
         }
         std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
         if (it == atom_indices_map.end()) {
            std::vector<int> v;
            atom_indices_map[alt_loc] = v;
         }
      }
   }

   int link_size = geom.link_size();
   for (int i = 0; i < link_size; i++) {
      if (geom.link(i).link_id == link_type) {
         unsigned int n_plane_restraints = geom.link(i).link_plane_restraint.size();
         for (unsigned int ip = 0; ip < n_plane_restraints; ip++) {
            std::vector<bool> fixed_flags(geom.link(i).link_plane_restraint[ip].n_atoms(), false);

            // Clear all atom index vectors
            for (auto& kv : atom_indices_map) {
               kv.second.clear();
            }

            unsigned int n_plane_atoms = geom.link(i).link_plane_restraint[ip].n_atoms();

            for (unsigned int irest_at = 0; irest_at < n_plane_atoms; irest_at++) {
               const auto& plane_restraint = geom.link(i).link_plane_restraint[ip];

               int link_res_n_atoms;
               const gemmi::Residue* res;

               if (plane_restraint.atom_comp_ids[irest_at] == 1) {
                  link_res_n_atoms = n_first_res_atoms;
                  res = first.residue;
                  fixed_flags[irest_at] = is_fixed_first_res;
               } else {
                  link_res_n_atoms = n_second_res_atoms;
                  res = second.residue;
                  fixed_flags[irest_at] = is_fixed_second_res;
               }

               if (!res) continue;

               std::string target_atom_name = plane_restraint.atom_id(irest_at);

               for (int iat = 0; iat < link_res_n_atoms; iat++) {
                  const auto& atom = res->atoms[iat];
                  // Pad atom name to 4-char MMDB format for comparison with atom_id()
                  std::string pdb_atom_name = coot::atom_id_mmdb_expand(atom.name);

                  if (target_atom_name == pdb_atom_name) {
                     std::string key;
                     if (atom.altloc != '\0') {
                        key = std::string(1, atom.altloc);
                     }
                     int idx_t_2 = iat;
                     atom_indices_map[key].push_back(idx_t_2);
                  }
               }
            }

            for (auto& kv : atom_indices_map) {
               if ((int)kv.second.size() > 3) {

                  // In the real code, make_fixed_flags returns flags based on
                  // bonded atoms — we replicate by defaulting to all false
                  std::vector<bool> other_fixed_flags(kv.second.size(), false);

                  for (unsigned int ii = 0; ii < other_fixed_flags.size(); ii++)
                     if (other_fixed_flags[ii])
                        fixed_flags[ii] = true;

                  // position and sigma
                  std::vector<std::pair<int, double> > pos_sigma;
                  for (unsigned int ii = 0; ii < kv.second.size(); ii++) {
                     double sigma_esd = geom.link(i).link_plane_restraint[ip].dist_esd();
                     pos_sigma.push_back(std::make_pair(kv.second[ii], sigma_esd));
                  }

                  gemmi_plane_restraint_data_t prd;
                  prd.pos_sigma = pos_sigma;
                  prd.fixed_flags = fixed_flags;
                  planes_out.push_back(prd);

                  n_plane++;
               }
            }
         }
      }
   }

   return n_plane;
}

} // namespace coot