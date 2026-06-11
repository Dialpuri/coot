#include "function.hh"
#include <gemmi/math.hpp>

#include "geometry/hb-types.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__h_bonds__mark_donors_and_acceptors/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__h_bonds__make_neighbour_map/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__h_bonds__make_h_bond_from_ligand_hydrogen/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__h_bonds__make_h_bond_from_environment_residue_hydrogen/gemmi/function.hh"

namespace coot {

// Helper: squared distance between two atoms
static double sq_dist(const gemmi::Atom& a, const gemmi::Atom& b) {
  double dx = a.pos.x - b.pos.x;
  double dy = a.pos.y - b.pos.y;
  double dz = a.pos.z - b.pos.z;
  return dx*dx + dy*dy + dz*dz;
}

// Convert ligand_hydrogen h_bond_gemmi (from dep port, in coot::h_bonds namespace) to our mcdonald_h_bond
static mcdonald_h_bond to_mcdonald(const h_bond_gemmi& src) {
  mcdonald_h_bond dst;
  dst.dist = src.dist;
  dst.angle_1 = src.angle_1;
  dst.angle_2 = src.angle_2;
  dst.angle_3 = src.angle_3;
  dst.hb_hydrogen = src.hb_hydrogen;
  dst.acceptor = src.acceptor;
  dst.donor = src.donor;
  dst.donor_neigh = src.donor_neigh;
  dst.acceptor_neigh = src.acceptor_neigh;
  dst.ligand_atom_is_donor = src.ligand_atom_is_donor;
  dst.bond_has_hydrogen_flag = (src.hb_hydrogen != nullptr);
  return dst;
}

// Convert environment_residue h_bond_gemmi (from dep port, in coot::h_bonds_gemmi namespace) to our mcdonald_h_bond
static mcdonald_h_bond to_mcdonald_env(const coot::h_bonds_gemmi::h_bond_gemmi& src) {
  mcdonald_h_bond dst;
  dst.dist = src.dist;
  dst.angle_1 = src.angle_1;
  dst.angle_2 = src.angle_2;
  dst.angle_3 = src.angle_3;
  dst.hb_hydrogen = src.hb_hydrogen;
  dst.acceptor = src.acceptor;
  dst.donor = src.donor;
  dst.donor_neigh = src.donor_neigh;
  dst.acceptor_neigh = src.acceptor_neigh;
  dst.ligand_atom_is_donor = src.ligand_atom_is_donor;
  dst.bond_has_hydrogen_flag = (src.hb_hydrogen != nullptr);
  return dst;
}

std::vector<mcdonald_h_bond>
get_mcdonald_and_thornton_gemmi(
    const std::vector<gemmi::CRA>& sel_1_atoms,
    const std::vector<gemmi::CRA>& sel_2_atoms,
    gemmi::Model& model,
    const coot::protein_geometry& geom,
    int imol,
    double max_dist) {

  bool debug = false;
  std::vector<mcdonald_h_bond> v;

  // mark_donors_and_acceptors_gemmi returns a map of atom* → hb_type
  std::map<gemmi::Atom*, int> hb_types_map =
    coot::h_bonds::mark_donors_and_acceptors_gemmi(sel_1_atoms, sel_2_atoms, geom, imol);

  // These distance are from the acceptor to the H - not the donor
  double min_dist = 0.1; // H-bonds are longer than this

  // Build neighbour map using existing _gemmi port
  std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>> neighbour_map =
    coot::h_bonds::make_neighbour_map_gemmi(sel_1_atoms, sel_2_atoms);

  // Build contact list: pairs within min_dist..max_dist
  struct Contact {
    int id1;
    int id2;
  };
  std::vector<Contact> contacts;
  double min_dist_sq = min_dist * min_dist;
  double max_dist_sq = max_dist * max_dist;

  for (int i = 0; i < static_cast<int>(sel_1_atoms.size()); i++) {
    if (!sel_1_atoms[i].atom) continue;
    for (int j = 0; j < static_cast<int>(sel_2_atoms.size()); j++) {
      if (!sel_2_atoms[j].atom) continue;
      double d_sq = sq_dist(*sel_1_atoms[i].atom, *sel_2_atoms[j].atom);
      if (d_sq >= min_dist_sq && d_sq <= max_dist_sq) {
        Contact c;
        c.id1 = i;
        c.id2 = j;
        contacts.push_back(c);
      }
    }
  }

  if (!contacts.empty()) {
    for (const auto& contact : contacts) {
      const gemmi::CRA& cra_1 = sel_1_atoms[contact.id1];
      const gemmi::CRA& cra_2 = sel_2_atoms[contact.id2];

      const gemmi::Atom* at_1 = cra_1.atom;
      const gemmi::Atom* at_2 = cra_2.atom;

      if (!at_1 || !at_2) continue;

      // Check alt_loc mismatch — skip if different non-empty alt locs
      char alt_1 = at_1->altloc;
      char alt_2 = at_2->altloc;
      if (alt_1 != '\0' && alt_2 != '\0' && alt_1 != alt_2)
        continue;

      // They must be on different residues
      if (cra_1.residue == cra_2.residue)
        continue;

      int hb_type_1 = coot::hb_t::HB_UNASSIGNED;
      int hb_type_2 = coot::hb_t::HB_UNASSIGNED;

      auto it1 = hb_types_map.find(const_cast<gemmi::Atom*>(at_1));
      if (it1 != hb_types_map.end()) hb_type_1 = it1->second;
      auto it2 = hb_types_map.find(const_cast<gemmi::Atom*>(at_2));
      if (it2 != hb_types_map.end()) hb_type_2 = it2->second;

      if (debug)
        std::cout << "DEBUG:: in get_mcdonald_and_thornton() "
                  << at_1->name << " " << at_2->name
                  << " hb_type: " << hb_type_1 << " " << hb_type_2 << std::endl;

      // hydrogen on ligand
      if (hb_type_1 == coot::hb_t::HB_HYDROGEN) {
        if (hb_type_2 == coot::hb_t::HB_ACCEPTOR || hb_type_2 == coot::hb_t::HB_BOTH) {
          // Get neighbour lists
          std::vector<std::pair<const gemmi::Atom*, float>> nb_1;
          std::vector<std::pair<const gemmi::Atom*, float>> nb_2;

          auto nit1 = neighbour_map.find(const_cast<gemmi::Atom*>(at_1));
          if (nit1 != neighbour_map.end()) {
            for (const auto& p : nit1->second)
              nb_1.push_back({const_cast<const gemmi::Atom*>(p.first), p.second});
          }

          auto nit2 = neighbour_map.find(const_cast<gemmi::Atom*>(at_2));
          if (nit2 != neighbour_map.end()) {
            for (const auto& p : nit2->second)
              nb_2.push_back({const_cast<const gemmi::Atom*>(p.first), p.second});
          }

          auto b_hbond = coot::h_bonds::make_h_bond_from_ligand_hydrogen_gemmi(at_1, at_2, nb_1, nb_2);
          if (b_hbond.first)
            v.push_back(to_mcdonald(b_hbond.second));
        }
      }

      // hydrogen on environment (protein) residue
      // Allow a special alternative case where the acceptor
      // is on the ligand and the donor is a water (because
      // waters may not (probably do not) have hydrogens.
      if (hb_type_1 == coot::hb_t::HB_ACCEPTOR || hb_type_1 == coot::hb_t::HB_BOTH) {
        std::string res_name_2 = cra_2.residue->name;
        if (hb_type_2 == coot::hb_t::HB_HYDROGEN || res_name_2 == "HOH") {
          std::vector<std::pair<const gemmi::Atom*, float>> nb_1;
          std::vector<std::pair<const gemmi::Atom*, float>> nb_2;

          auto nit1 = neighbour_map.find(const_cast<gemmi::Atom*>(at_1));
          if (nit1 != neighbour_map.end()) {
            for (const auto& p : nit1->second)
              nb_1.push_back({const_cast<const gemmi::Atom*>(p.first), p.second});
          }

          auto nit2 = neighbour_map.find(const_cast<gemmi::Atom*>(at_2));
          if (nit2 != neighbour_map.end()) {
            for (const auto& p : nit2->second)
              nb_2.push_back({const_cast<const gemmi::Atom*>(p.first), p.second});
          }

          auto b_hbond = coot::h_bonds_gemmi::make_h_bond_from_environment_residue_hydrogen_gemmi(at_1, at_2, nb_1, nb_2);
          if (b_hbond.first) {
            if (debug)
              std::cout << "DEBUG:: ===> in get_m&d: pushing back b_hbond\n" << std::endl;
            v.push_back(to_mcdonald_env(b_hbond.second));
          } else {
            if (debug)
              std::cout << "DEBUG:: reject\n" << std::endl;
          }
        }
      }
    }
  }

  return v;
}

} // namespace coot