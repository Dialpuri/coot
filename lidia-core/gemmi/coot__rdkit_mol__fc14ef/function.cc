#include "function.hh"

#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "lidia-core/rdkit-interface.hh"
#include "utils/coot-utils.hh"
#include "lidia-core/lbg-molfile.hh"

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <iomanip>

// Helper: pad atom name to 4 chars in PDB/CCP4/MMDB format for comparison with restraints
static std::string pad4(const std::string &s) {
   std::string trimmed = s;
   size_t a = trimmed.find_first_not_of(' ');
   size_t b = trimmed.find_last_not_of(' ');
   if (a == std::string::npos) return "    ";
   trimmed = trimmed.substr(a, b - a + 1);
   std::string p = ' ' + trimmed;
   while (p.size() < 4) p += ' ';
   return p;
}

// Gemmi-adapted version of get_chiral_tag (from rdkit-interface.cc)
static RDKit::Atom::ChiralType
get_chiral_tag_gemmi(gemmi::Residue *residue_p,
                     const coot::dictionary_residue_restraints_t &restraints,
                     gemmi::Atom *atom_p) {
   RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED;
   if (!residue_p) return chiral_tag;
   std::string atom_name = pad4(atom_p->name);

   for (unsigned int ichi=0; ichi<restraints.chiral_restraint.size(); ichi++) {
      if (restraints.chiral_restraint[ichi].atom_id_c_4c() == atom_name) {
         const coot::dict_chiral_restraint_t &chiral_restraint = restraints.chiral_restraint[ichi];
         int n_neigbours_found = 0;
         std::vector<int> ni(4, -1);
         for (unsigned int iat=0; iat<residue_p->atoms.size(); iat++) {
            std::string an = pad4(residue_p->atoms[iat].name);
            if (an == chiral_restraint.atom_id_1_4c()) { ni[1] = static_cast<int>(iat); n_neigbours_found++; }
            if (an == chiral_restraint.atom_id_2_4c()) { ni[2] = static_cast<int>(iat); n_neigbours_found++; }
            if (an == chiral_restraint.atom_id_3_4c()) { ni[3] = static_cast<int>(iat); n_neigbours_found++; }
         }
         if (n_neigbours_found == 3) {
            bool atom_orders_match = false;
            if ((ni[3] > ni[2]) && (ni[2] > ni[1])) atom_orders_match = true;
            if ((ni[1] > ni[3]) && (ni[3] > ni[2])) atom_orders_match = true;
            if ((ni[2] > ni[1]) && (ni[1] > ni[3])) atom_orders_match = true;
            if (atom_orders_match) {
               if (chiral_restraint.volume_sign == 1) chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
               if (chiral_restraint.volume_sign == -1) chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
            } else {
               if (chiral_restraint.volume_sign == 1) chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
               if (chiral_restraint.volume_sign == -1) chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
            }
         }
         break;
      }
   }
   return chiral_tag;
}

// Gemmi-adapted set_atom_chirality
static void set_atom_chirality_gemmi(RDKit::Atom *rdkit_at,
                                      gemmi::Atom *atom,
                                      gemmi::Residue *residue_p,
                                      const coot::dictionary_residue_restraints_t &restraints) {
    std::string an = pad4(atom->name);
    for (unsigned int i = 0; i < restraints.atom_info.size(); i++) {
        if (restraints.atom_info[i].atom_id_4c == an) {
            RDKit::Atom::ChiralType chiral_tag = get_chiral_tag_gemmi(residue_p, restraints, atom);
            if (chiral_tag != RDKit::Atom::CHI_UNSPECIFIED) {
                rdkit_at->setChiralTag(chiral_tag);
            }
            break;
        }
    }
}

RDKit::RWMol coot::rdkit_mol_gemmi(
    gemmi::CRA cra,
    const coot::dictionary_residue_restraints_t &restraints,
    const std::string &alt_conf,
    bool do_undelocalize) {

   if (cra.residue == nullptr) {
      throw std::runtime_error("rdkit_mol_gemmi: null residue");
   }
   if (cra.chain == nullptr) {
      throw std::runtime_error("rdkit_mol_gemmi: null chain");
   }

   auto get_alt_confs_in_residue = [](gemmi::Residue *residue_p) {
      std::vector<std::string> v;
      for (auto &atom : residue_p->atoms) {
         std::string ac(1, atom.altloc);
         if (std::find(v.begin(), v.end(), ac) == v.end()) {
            v.push_back(ac);
         }
      }
      return v;
   };

   bool debug = false;

   if (debug) {
      coot::residue_spec_t spec(cra.chain->name, cra.residue->seqid.num.value,
                                cra.residue->seqid.icode != ' ' ? std::string(1, cra.residue->seqid.icode) : "");
      std::cout << "=========== in rdkit_mol_gemmi() with restraints that have "
                << restraints.atom_info.size() << " atoms, "
                << restraints.bond_restraint.size() << " bond restraints with do_undelocalize "
                << do_undelocalize
                << " for residue " << spec.format()
                << " and alt conf " << "\"" << alt_conf << "\"" << std::endl;
   }

   if (debug)
      for (unsigned int ii=0; ii<restraints.atom_info.size(); ii++)
         std::cout << ii << "   " << restraints.atom_info[ii] << std::endl;

   RDKit::RWMol m;

   std::string n = coot::util::remove_trailing_whitespace(restraints.residue_info.name);
   m.setProp("_Name", n);
   m.setProp("ResName", cra.residue->name);
   m.setProp("ResNumber", cra.residue->seqid.num.value);
   m.setProp("ChainID", cra.chain->name);
   m.setProp("alt_id", alt_conf);

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   
   std::vector<gemmi::Atom *> residue_atoms;
   for (auto &at : cra.residue->atoms) {
      std::string padded_name = pad4(at.name);
      if (padded_name == " OXT" || padded_name == " HN ")
         continue;
      residue_atoms.push_back(&at);
   }
   int n_residue_atoms = static_cast<int>(residue_atoms.size());

   std::vector<std::string> added_atom_names;
   std::vector<gemmi::Atom *> added_atoms;
   std::map<std::string, int> atom_index;
   int current_atom_id = 0;
   std::vector<std::pair<int, int> > bonded_atoms;
  
   for (int iat_1=0; iat_1<n_residue_atoms; iat_1++) {
      gemmi::Atom *at_1 = residue_atoms[iat_1];
      std::string atom_name_1 = pad4(at_1->name);
      
      bool found_a_bonded_atom = false;
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
         if (restraints.bond_restraint[ib].atom_id_1_4c() == atom_name_1) {
            for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
               if (pad4(residue_atoms[iat_2]->name) == restraints.bond_restraint[ib].atom_id_2_4c()) {
                  found_a_bonded_atom = true;
                  break;
               } 
            }
         } 
         if (restraints.bond_restraint[ib].atom_id_2_4c() == atom_name_1) {
            for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
               if (pad4(residue_atoms[iat_2]->name) == restraints.bond_restraint[ib].atom_id_1_4c()) {
                  found_a_bonded_atom = true;
                  break;
               }
            }
         }
      }

      int ai_idx = -1;
      for (unsigned int jj=0; jj<restraints.atom_info.size(); jj++) {
         if (restraints.atom_info[jj].atom_id_4c == atom_name_1) {
            ai_idx = jj;
            break;
         }
      }

      if (found_a_bonded_atom) {
         if (ai_idx > -1) {
            std::pair<int, int> p(iat_1, ai_idx);
            bool already_there = false;
            for (unsigned int iba=0; iba<bonded_atoms.size(); iba++) {
               if (bonded_atoms[iba].second == ai_idx) {
                  already_there = true;
                  break;
               }
            }
            if (! already_there) {
               bonded_atoms.push_back(p);
            } else {
               std::cout << "WARNING:: rdkit_mol_gemmi() atom index already there " << ai_idx << std::endl;
            }
         }
      }
   }

   if (debug) { 
      std::cout << "DEBUG:: number of bonded atoms with alt conf \"" << alt_conf << "\" found: "
                << bonded_atoms.size() << std::endl;
   }

   if (! bonded_atoms.empty()) {

      for (unsigned int iat=0; iat<bonded_atoms.size(); iat++) {
         gemmi::Atom *at = residue_atoms[bonded_atoms[iat].first];
         std::string atom_name = pad4(at->name);
         
         if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name) != added_atom_names.end()) {
            std::cout << "!!!! Problem? atom name \"" << atom_name
                      << "\" was already added" << std::endl;
         } else {
            RDKit::Atom *rdkit_at = new RDKit::Atom;
            try {
               std::string ele_capped =
                  coot::util::capitalise(coot::util::remove_leading_spaces(at->element.name()));
               int atomic_number = tbl->getAtomicNumber(ele_capped);
               rdkit_at->setAtomicNum(atomic_number);
               rdkit_at->setIsotope(0);
               rdkit_at->setProp("name", atom_name);

               const coot::dict_atom &atom_info = restraints.atom_info[bonded_atoms[iat].second];
               if (atom_info.formal_charge.first)
                  rdkit_at->setFormalCharge(atom_info.formal_charge.second);

               std::string type_energy = restraints.type_energy(atom_name);
               if (type_energy != "") {
                  if (type_energy == "NT") {
                     bool charge_it = true;
                     bool include_H_neighb_bonds = false;
                     if (restraints.neighbours(atom_name, include_H_neighb_bonds).size() == 3)
                        charge_it = false;
                     if (charge_it)
                        rdkit_at->setFormalCharge(1);
                  }
               }

               set_atom_chirality_gemmi(rdkit_at, at, cra.residue, restraints);
               m.addAtom(rdkit_at);
               
               added_atom_names.push_back(atom_name);
               added_atoms.push_back(residue_atoms[iat]);
               atom_index[atom_name] = current_atom_id;
               current_atom_id++;
            }
            catch (const std::exception &rte) {
               std::cout << rte.what() << std::endl;
            }
         }
      }

      if (debug) {
         std::cout << "DEBUG:: number of atoms in rdkit mol: " << m.getNumAtoms() << std::endl;
      } 

      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
         RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
         RDKit::Bond *bond = new RDKit::Bond(type);

         if (restraints.bond_restraint[ib].type() == "deloc") {
            bond->setProp("restraints-type", "deloc");
         }
      
         std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
         std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
         std::string ele_1 = restraints.element(atom_name_1);
         std::string ele_2 = restraints.element(atom_name_2);
         int idx_1 = -1;
         int idx_2 = -1;

         for (unsigned int iat=0; iat<added_atom_names.size(); iat++) {
            if (added_atom_names[iat] == atom_name_1) { idx_1 = iat; break; }
         }
         for (unsigned int iat=0; iat<added_atom_names.size(); iat++) { 
            if (added_atom_names[iat] == atom_name_2) { idx_2 = iat; break; }
         }

         if (idx_1 != -1) { 
            if (idx_2 != -1) {
               bool swap_order = false;
               if (restraints.chiral_restraint.size()) {
                  swap_order = chiral_check_order_swap(m[idx_1], m[idx_2], restraints.chiral_restraint);
               } else {
                  swap_order = chiral_check_order_swap(m[idx_1], m[idx_2]);
               }
               if (! swap_order) {
                  swap_order = chiral_check_order_swap_singleton(m[idx_1], m[idx_2], restraints);
               }
               if (! swap_order) {
                  bond->setBeginAtomIdx(idx_1);
                  bond->setEndAtomIdx(idx_2);
               } else {
                  bond->setBeginAtomIdx(idx_2);
                  bond->setEndAtomIdx(idx_1);
               } 
               if (type == RDKit::Bond::AROMATIC) {
                  bond->setIsAromatic(true);
                  m[idx_1]->setIsAromatic(true);
                  m[idx_2]->setIsAromatic(true);
               }
               m.addBond(bond);
            } else {
               if (ele_2 != " H") {
                  if (atom_name_2 != " OXT" && atom_name_1 != " O1 ") {
                     std::cout << "WARNING:: oops, bonding in rdkit_mol_gemmi() "
                               << "failed to get atom index idx_2 for atom name: "
                               << atom_name_2 << " ele :" << ele_2 << ":" << std::endl;
                  }
                  std::string message = "WARNING:: rdkit_mol_gemmi(): Failed to get atom index for atom name \"";
                  message += atom_name_2;
                  message += "\" in residue of type " + cra.residue->name;
                  throw std::runtime_error(message);
               }
            }
         } else {
            if (ele_1 != " H") { 
               if (atom_name_2 != " OXT" && atom_name_1 != " O1 ") {
                  std::cout << "WARNING:: oops, bonding in rdkit_mol_gemmi() "
                            << "failed to get atom index idx_1 for atom name: \""
                            << atom_name_1 << "\" ele :" << ele_1 << ":" << std::endl;
               }
               std::string message = "WARNING:: rdkit_mol_gemmi(): Failed to get atom index for atom name \"";
               message += atom_name_1;
               message += "\" in residue of type " + cra.residue->name;
               throw std::runtime_error(message);
            }
         }
         delete bond;
      }

      if (debug) { 
         std::cout << "DEBUG:: rdkit_mol_gemmi() number of bond restraints:    "
                   << restraints.bond_restraint.size() << std::endl;
         debug_rdkit_molecule(&m);
      }

      std::vector<int> Hs_added_list;
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 
         RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
         if (type == RDKit::Bond::AROMATIC) {
            std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
            std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
            int idx_1 = -1;
            int idx_2 = -1;
            for (unsigned int iat=0; iat<m.getNumAtoms(); iat++) {
               try {
                  std::string name;
                  RDKit::Atom *at_p = m[iat];
                  at_p->getProp("name", name);
                  if (name == atom_name_1) idx_1 = iat;
                  if (name == atom_name_2) idx_2 = iat;
               }
               catch (const KeyErrorException &) {}
            }
            if (idx_1 != -1 && idx_2 != -1) {         
               if (m[idx_1]->getAtomicNum() == 7) {
                  if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_1) == Hs_added_list.end()) { 
                     std::string n = add_H_to_ring_N_as_needed(&m, idx_1, atom_name_1, restraints);
                     if (n != "") added_atom_names.push_back(n);
                     Hs_added_list.push_back(idx_1);
                  }
               }
               if (m[idx_2]->getAtomicNum() == 7) {
                  if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_2) == Hs_added_list.end()) { 
                     std::string n = add_H_to_ring_N_as_needed(&m, idx_2, atom_name_2, restraints);
                     if (n != "") added_atom_names.push_back(n);
                     Hs_added_list.push_back(idx_2);
                  }
               }
            }
         }
      }

      if (debug) debug_rdkit_molecule(&m);
      if (do_undelocalize) coot::undelocalise(&m);
      coot::assign_formal_charges(&m);

      std::vector<std::vector<int> > ring_info;
      RDKit::MolOps::findSSSR(m, ring_info);
      RDKit::MolOps::cleanUp(m);
      RDKit::MolOps::sanitizeMol(m);

      unsigned int n_atoms = m.getNumAtoms();
      if (n_atoms > 0) {
         RDKit::UINT_VECT ranks(m.getNumAtoms(), -1);
         RDKit::Chirality::assignAtomCIPRanks(m, ranks);
         for (unsigned int iat=0; iat<bonded_atoms.size(); iat++) {
            const coot::dict_atom &atom_info = restraints.atom_info[bonded_atoms[iat].second];
            if (atom_info.pdbx_stereo_config.first &&
                (atom_info.pdbx_stereo_config.second == "R" ||
                 atom_info.pdbx_stereo_config.second == "S")) {
               std::vector<std::pair<const RDKit::Atom *, unsigned int> > neighbs;
               unsigned int idx_iat = bonded_atoms[iat].first;
               if (idx_iat >= n_atoms) {
                  std::cout << "ERROR:: rdkit_mol_gemmi() chiral-check: trying to get atom with "
                            << "index " << idx_iat << " but molecule has " << n_atoms << " atoms" << std::endl;
               } else {
                  RDKit::Atom *rdkit_at = m[idx_iat];
                  RDKit::ROMol::OEDGE_ITER beg,end;
                  boost::tie(beg,end) = m.getAtomBonds(rdkit_at);
                  while(beg != end){
                     const RDKit::Bond *bond=m[*beg]; ++beg;
                     const RDKit::Atom *nbr=bond->getOtherAtom(rdkit_at);
                     unsigned int cip_rank = 0;
                     nbr->getProp(RDKit::common_properties::_CIPRank, cip_rank);
                     neighbs.push_back(std::pair<const RDKit::Atom *, unsigned int>(nbr, cip_rank));
                  }
                  std::vector<std::pair<const RDKit::Atom *, unsigned int> > sorted_neighbs = neighbs;
                  std::sort(sorted_neighbs.begin(), sorted_neighbs.end(), cip_rank_sorter);
                  bool inverted = true;
                  if (neighbs.size() == 3) {
                     neighbs.resize(3); sorted_neighbs.resize(3);
                     if (neighbs[0]==sorted_neighbs[0] && neighbs[1]==sorted_neighbs[1] && neighbs[2]==sorted_neighbs[2]) inverted = false;
                     if (neighbs[0]==sorted_neighbs[1] && neighbs[1]==sorted_neighbs[2] && neighbs[2]==sorted_neighbs[0]) inverted = false;
                     if (neighbs[0]==sorted_neighbs[2] && neighbs[1]==sorted_neighbs[0] && neighbs[2]==sorted_neighbs[1]) inverted = false;
                  } else if (neighbs.size() == 4) {
                     std::vector<const RDKit::Atom *> needed_atoms(3);
                     needed_atoms[0] = sorted_neighbs[1].first;
                     needed_atoms[1] = sorted_neighbs[2].first;
                     needed_atoms[2] = sorted_neighbs[3].first;
                     unsigned int n_found = 0;
                     for (unsigned int jj=0; jj<3; jj++)
                        for (unsigned int ii=0; ii<3; ii++)
                           if (needed_atoms[ii] == neighbs[jj].first) n_found += 1;
                     if (n_found == 3) {
                        if (neighbs[0]==sorted_neighbs[0] && neighbs[1]==sorted_neighbs[1] && neighbs[2]==sorted_neighbs[2]) inverted = false;
                        if (neighbs[0]==sorted_neighbs[1] && neighbs[1]==sorted_neighbs[2] && neighbs[2]==sorted_neighbs[0]) inverted = false;
                        if (neighbs[0]==sorted_neighbs[2] && neighbs[1]==sorted_neighbs[0] && neighbs[2]==sorted_neighbs[1]) inverted = false;
                        if (atom_info.pdbx_stereo_config.second == "R") {
                           if (inverted) rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                           else rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                        }
                        if (atom_info.pdbx_stereo_config.second == "S") {
                           if (inverted) rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                           else rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                        }
                     } else {
                        unsigned int idx_cip_rank_lowest = 0;
                        unsigned int cip_rank_lowest = 99999;
                        for (unsigned int jj=0; jj<4; jj++) {
                           if (neighbs[jj].second < cip_rank_lowest) {
                              cip_rank_lowest = neighbs[jj].second;
                              idx_cip_rank_lowest = jj;
                           }
                        }
                        if (idx_cip_rank_lowest == 1) inverted = true;
                        if (idx_cip_rank_lowest == 2) inverted = false;
                        if (idx_cip_rank_lowest == 3) inverted = true;
                     }
                  } else {
                     std::cout << "WARNING:: crazy atom - too many connections " << atom_info << std::endl;
                  }
                  if (atom_info.pdbx_stereo_config.second == "R") {
                     if (inverted) rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                     else rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                  }
                  if (atom_info.pdbx_stereo_config.second == "S") {
                     if (inverted) rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
                     else rdkit_at->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
                  }
               }
            }
         }

         std::vector<std::string> alt_confs_in_residue = get_alt_confs_in_residue(cra.residue);
         for (unsigned int iconf=0; iconf<alt_confs_in_residue.size(); iconf++) {
            RDKit::Conformer *conf = new RDKit::Conformer(m.getNumAtoms());
            conf->set3D(true);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               std::string atom_name = pad4(residue_atoms[iat]->name);
               std::map<std::string, int>::const_iterator it = atom_index.find(atom_name);
               if (it != atom_index.end()) {
                  RDGeom::Point3D pos(residue_atoms[iat]->pos.x,
                                      residue_atoms[iat]->pos.y,
                                      residue_atoms[iat]->pos.z);
                  conf->setAtomPos(it->second, pos);
               }
            }
            int conf_id = m.addConformer(conf);
            RDKit::MolOps::assignChiralTypesFrom3D(m, conf_id, true);
         }
         set_energy_lib_atom_types(&m);
      }
   }
   
   return m;
}