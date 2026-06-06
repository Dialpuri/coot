#include "function.hh"
#include <gemmi/neighbor.hpp>
#include <gemmi/calculate.hpp>
#include <cmath>
#include <algorithm>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__contains_nucleic_acid/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__mark_donors_and_acceptors/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__is_nucleotide/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__similar_normals/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

// Pad atom name (same as in mark_donors_and_acceptors)
static inline std::string pad_atom_name(const std::string& name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    if (name.size() == 3) return " " + name;
    return name;
}

std::vector<coot::stack_and_pair::paired_residues_info_t>
coot::stack_and_pair::paired_residues_gemmi(
    gemmi::Model &model,
    const std::vector<std::pair<bool, gemmi::Residue *>> &residues_vec,
    bool residues_are_all_moving_flag,
    const coot::protein_geometry &geom,
    int imol_enc) {

   std::vector<coot::stack_and_pair::paired_residues_info_t> v;

   float dist_crit = 3.2;
   float dist_crit_sqrt = dist_crit * dist_crit;

   // Collect all atom CRAs from the model (selected_atoms_all equivalent)
   std::vector<gemmi::CRA> all_atom_cras;
   for (auto &chain : model.chains) {
      for (auto &res : chain.residues) {
         for (auto &atom : res.atoms) {
            all_atom_cras.push_back(gemmi::CRA{&chain, &res, &atom});
         }
      }
   }

   // Collect moving atom CRAs
   std::vector<gemmi::CRA> moving_atom_cras;
   // Hack from original: residues_are_all_moving_flag is forced to true
   residues_are_all_moving_flag = true;

   if (residues_are_all_moving_flag) {
      moving_atom_cras = all_atom_cras;
   } else {
      // When not all moving, collect atoms from the non-moving residues only
      for (const auto &pr : residues_vec) {
         if (!pr.first) continue;
         if (!pr.second) continue;
         for (auto &atom : pr.second->atoms) {
            // We need the chain pointer - find it
            for (auto &chain : model.chains) {
               for (auto &res : chain.residues) {
                  if (&res == pr.second) {
                     moving_atom_cras.push_back(gemmi::CRA{&chain, &res, &atom});
                  }
               }
            }
         }
      }
   }

   // Build a vector of moving residues for contains_nucleic_acid check
   std::vector<gemmi::Residue> moving_residues;
   std::set<gemmi::Residue*> seen_residues;
   for (const auto &cra : moving_atom_cras) {
      if (cra.residue && seen_residues.find(cra.residue) == seen_residues.end()) {
         seen_residues.insert(cra.residue);
         moving_residues.push_back(*cra.residue);
      }
   }

   if (static_cast<int>(moving_atom_cras.size()) > 1) {
      // Check for nucleic acid
      if (coot::stack_and_pair::contains_nucleic_acid_gemmi(moving_residues)) {

         // Build normal_map from the residues
         // For now empty - similar_normals requires it
         std::map<gemmi::Residue *, clipper::Coord_orth> normal_map;

         // Mark donors and acceptors for ALL atoms
         std::vector<int> hb_types =
             coot::stack_and_pair::mark_donors_and_acceptors_gemmi(all_atom_cras, geom, imol_enc);

         // Build a lookup: Atom pointer -> index in all_atom_cras
         std::map<gemmi::Atom const*, size_t> atom_to_idx;
         for (size_t i = 0; i < all_atom_cras.size(); i++) {
            if (all_atom_cras[i].atom)
               atom_to_idx[all_atom_cras[i].atom] = i;
         }

         // Use gemmi::NeighborSearch to find contacts within dist_crit
         gemmi::NeighborSearch ns(model, gemmi::UnitCell(), static_cast<double>(dist_crit));
         ns.populate(false); // include_h = false

         // interesting Hydrogen bonds - but not base pairing
         std::set<std::string> excluded_oxygens;
         excluded_oxygens.insert(" OP1"); excluded_oxygens.insert(" OP2");
         excluded_oxygens.insert(" O2'"); excluded_oxygens.insert(" O3'");
         excluded_oxygens.insert(" O5'");

         std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_vec;

         // For each moving atom, find nearby atoms in all selection
         for (const auto &m_cra : moving_atom_cras) {
            if (!m_cra.atom) continue;
            gemmi::Atom *at_1 = m_cra.atom;

            // Check element is O or N (gemmi: element.name())
            std::string ele_1 = at_1->element.name();
            if (ele_1 != "O" && ele_1 != "N") continue;

            auto hits = ns.find_atoms(at_1->pos, '\0', 0.0, static_cast<double>(dist_crit));

            for (auto *m : hits) {
               gemmi::CRA hit_cra = m->to_cra(model);
               if (!hit_cra.atom) continue;
               gemmi::Atom *at_2 = hit_cra.atom;

               // in different residues
               if (m_cra.residue != hit_cra.residue) {
                  std::string ele_2 = at_2->element.name();
                  if (ele_2 != "O" && ele_2 != "N") continue;

                  // distance check (already done by find_atoms with radius, but verify)
                  double dx = at_1->pos.x - at_2->pos.x;
                  double dy = at_1->pos.y - at_2->pos.y;
                  double dz = at_1->pos.z - at_2->pos.z;
                  double dd = dx * dx + dy * dy + dz * dz;
                  if (dd >= dist_crit_sqrt) continue;

                  // Get hb types from the parallel vector
                  int hb_type_1 = 0;
                  int hb_type_2 = 0;

                  auto it1 = atom_to_idx.find(at_1);
                  auto it2 = atom_to_idx.find(at_2);

                  if (it1 != atom_to_idx.end() && it1->second < hb_types.size()) {
                     hb_type_1 = hb_types[it1->second];
                  }
                  if (it2 != atom_to_idx.end() && it2->second < hb_types.size()) {
                     hb_type_2 = hb_types[it2->second];
                  }

                  // Donor-acceptor check
                  if ((hb_type_1 == coot::hb_t::HB_ACCEPTOR || hb_type_1 == coot::hb_t::HB_BOTH)) {
                     if ((hb_type_2 == coot::hb_t::HB_DONOR || hb_type_2 == coot::hb_t::HB_BOTH)) {

                        // Same chain check - skip if too close in sequence on same chain
                        if (m_cra.chain && hit_cra.chain && m_cra.chain->name == hit_cra.chain->name) {
                           // Check residue index delta
                           int ri1 = 0, ri2 = 0;
                           bool found1 = false, found2 = false;
                           {
                              for (const auto &c : model.chains) {
                                 if (c.name == m_cra.chain->name) {
                                    int ci = 0;
                                    for (const auto &r : c.residues) {
                                       if (&r == m_cra.residue) { ri1 = ci; found1 = true; break; }
                                       ci++;
                                    }
                                 }
                              }
                              for (const auto &c : model.chains) {
                                 if (c.name == hit_cra.chain->name) {
                                    int ci = 0;
                                    for (const auto &r : c.residues) {
                                       if (&r == hit_cra.residue) { ri2 = ci; found2 = true; break; }
                                       ci++;
                                    }
                                 }
                              }
                              if (found1 && found2 && std::abs(ri2 - ri1) < 2)
                                 continue;
                           }

                           // Both nucleotides?
                           if (coot::util::is_nucleotide_gemmi(*(m_cra.residue)) &&
                               coot::util::is_nucleotide_gemmi(*(hit_cra.residue))) {

                              // Similar normals check
                              if (coot::stack_and_pair::similar_normals_gemmi(
                                      m_cra.residue, hit_cra.residue, normal_map, 0)) {

                                 // Dot product checks
                                 clipper::Coord_orth pt_1 = coot::co_gemmi(at_1);
                                 clipper::Coord_orth pt_2 = coot::co_gemmi(at_2);
                                 clipper::Coord_orth diff(pt_2 - pt_1);
                                 clipper::Coord_orth atom_atom_unit_vector(diff.unit());

                                 double dp_1 = clipper::Coord_orth::dot(atom_atom_unit_vector, normal_map[m_cra.residue]);
                                 double dp_2 = clipper::Coord_orth::dot(atom_atom_unit_vector, normal_map[hit_cra.residue]);

                                 if (std::abs(dp_1) < 0.5 && std::abs(dp_2) < 0.5) {
                                    // Check excluded oxygens
                                    std::string name_1 = pad_atom_name(at_1->name);
                                    std::string name_2 = pad_atom_name(at_2->name);
                                    if (excluded_oxygens.find(name_1) == excluded_oxygens.end() &&
                                        excluded_oxygens.find(name_2) == excluded_oxygens.end()) {
                                       atom_vec.push_back(std::make_pair(at_1, at_2));
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

         // Convert atom_vec into paired_residues_info_t
         for (std::size_t i = 0; i < atom_vec.size(); i++) {
            gemmi::Atom *at_1 = atom_vec[i].first;
            gemmi::Atom *at_2 = atom_vec[i].second;
            gemmi::Residue *res_1 = nullptr;
            gemmi::Residue *res_2 = nullptr;

            // Find residues from CRA
            for (const auto &cra : moving_atom_cras) {
               if (cra.atom == at_1 && cra.residue) res_1 = cra.residue;
               if (cra.atom == at_2 && cra.residue) res_2 = cra.residue;
            }
            if (!res_1 || !res_2) continue;

            bool found = false;
            for (std::size_t j = 0; j < v.size(); j++) {
               if (v[j].res_1 == res_1 && v[j].res_2 == res_2) {
                  v[j].atom_pair_vec.push_back(atom_vec[i]);
                  found = true;
                  break;
               }
               if (v[j].res_1 == res_2 && v[j].res_2 == res_1) {
                  std::swap(atom_vec[i].first, atom_vec[i].second);
                  v[j].atom_pair_vec.push_back(atom_vec[i]);
                  found = true;
                  break;
               }
            }
            if (!found) {
               std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> apv;
               apv.push_back(std::make_pair(at_1, at_2));
               coot::stack_and_pair::paired_residues_info_t pri(res_1, res_2, apv);
               v.push_back(pri);
            }
         }
      }
   }

   return v;
}