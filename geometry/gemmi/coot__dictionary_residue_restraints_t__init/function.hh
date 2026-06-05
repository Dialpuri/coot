#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <clipper/core/coords.h>
#include <clipper/core/clipper_util.h>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cmath>
#include <algorithm>

namespace coot {

struct dict_chem_comp_t {
  std::string comp_id;
  std::string comp_name;
  std::string comp_full_name;
  std::string group;
  int n_atoms;
  int n_atoms_no_h;
  std::string description_level;

  dict_chem_comp_t() : n_atoms(0), n_atoms_no_h(0) {}
  dict_chem_comp_t(const std::string& a_comp_id,
                   const std::string& a_comp_name,
                   const std::string& a_comp_full_name,
                   std::string a_group,
                   int a_n_atoms,
                   int a_n_atoms_no_h,
                   const std::string& a_description_level)
    : comp_id(a_comp_id), comp_name(a_comp_name), comp_full_name(a_comp_full_name),
      group(a_group), n_atoms(a_n_atoms), n_atoms_no_h(a_n_atoms_no_h),
      description_level(a_description_level) {}
};

struct dict_atom {
  std::string atom_id;
  std::string altloc;
  std::string element;
  std::string name;
  std::pair<bool, clipper::Coord_orth> model_Cartn;
  std::pair<bool, float> model_Occ;
  std::pair<bool, float> model_B;

  dict_atom() : model_Cartn(std::make_pair(false, clipper::Coord_orth())),
                model_Occ(std::make_pair(false, 0)), model_B(std::make_pair(false, 0)) {}
  dict_atom(const std::string& a_atom_id, const std::string& a_altloc,
            const std::string& a_element, const std::string& a_name,
            std::pair<bool, clipper::Coord_orth> a_model_Cartn)
    : atom_id(a_atom_id), altloc(a_altloc), element(a_element), name(a_name),
      model_Cartn(a_model_Cartn), model_Occ(std::make_pair(false, 0)), model_B(std::make_pair(false, 0)) {}
};

struct dict_bond_restraint_t {
  std::string atom_id_1_;
  std::string atom_id_2_;
  std::string type_;
  double dist_;
  double dist_esd_;
  unsigned order_;
  bool is_hydrogen_bond_;
  double pdb_dist_;
  bool is_aromatic_;

  dict_bond_restraint_t() : dist_(0), dist_esd_(0), order_(0), is_hydrogen_bond_(false), pdb_dist_(0), is_aromatic_(false) {}

  dict_bond_restraint_t(std::string a_atom_id_1,
                        std::string a_atom_id_2,
                        std::string a_type,
                        double a_dist,
                        double a_dist_esd,
                        unsigned a_order,
                        bool a_is_hydrogen_bond,
                        double a_pdb_dist,
                        bool a_is_aromatic)
    : atom_id_1_(a_atom_id_1), atom_id_2_(a_atom_id_2), type_(a_type),
      dist_(a_dist), dist_esd_(a_dist_esd), order_(a_order),
      is_hydrogen_bond_(a_is_hydrogen_bond), pdb_dist_(a_pdb_dist), is_aromatic_(a_is_aromatic) {}

  std::string atom_id_1() const { return atom_id_1_; }
  std::string atom_id_2() const { return atom_id_2_; }
  std::string type() const { return type_; }
  double dist() const { return dist_; }
  double dist_esd() const { return dist_esd_; }
  unsigned order() const { return order_; }
  bool is_hydrogen_bond() const { return is_hydrogen_bond_; }
  double pdb_dist() const { return pdb_dist_; }
  bool is_aromatic() const { return is_aromatic_; }
};

struct dict_angle_restraint_t {
  std::string atom_id_1_;
  std::string atom_id_2_;
  std::string atom_id_3_;
  double angle_;
  double angle_esd_;

  dict_angle_restraint_t() : angle_(0), angle_esd_(0) {}

  dict_angle_restraint_t(std::string a_atom_id_1,
                         std::string a_atom_id_2,
                         std::string a_atom_id_3,
                         double a_angle,
                         double a_angle_esd)
    : atom_id_1_(a_atom_id_1), atom_id_2_(a_atom_id_2), atom_id_3_(a_atom_id_3),
      angle_(a_angle), angle_esd_(a_angle_esd) {}

  std::string atom_id_1() const { return atom_id_1_; }
  std::string atom_id_2() const { return atom_id_2_; }
  std::string atom_id_3() const { return atom_id_3_; }
  double angle() const { return angle_; }
  double angle_esd() const { return angle_esd_; }
};

class dictionary_residue_restraints_t {
public:
  std::vector<dict_bond_restraint_t> bond_restraint;
  std::vector<dict_angle_restraint_t> angle_restraint;
  std::vector<dict_atom> atom_info;
  dict_chem_comp_t residue_info;
  bool filled_with_bond_order_data_only_flag;
  bool nuclear_distances_flag;

  struct atom_pair_t {
    const gemmi::Atom *at_1;
    const gemmi::Atom *at_2;

    atom_pair_t() : at_1(nullptr), at_2(nullptr) {}
    atom_pair_t(const gemmi::Atom *a1, const gemmi::Atom *a2) : at_1(a1), at_2(a2) {}

    const gemmi::Atom *shared_atom(const atom_pair_t &other) const {
      if (at_1 == other.at_1) return at_1;
      if (at_1 == other.at_2) return at_1;
      if (at_2 == other.at_1) return at_2;
      if (at_2 == other.at_2) return at_2;
      return nullptr;
    }
  };

  dictionary_residue_restraints_t()
    : filled_with_bond_order_data_only_flag(false), nuclear_distances_flag(false) {}

  // Pad atom name to 4 characters in PDB format: " CA ", " N  ", etc.
  static std::string pad_name(const std::string& name) {
    std::string no_space;
    for (char c : name) { if (c != ' ') no_space += c; }
    return " " + no_space + std::string(std::max(0, 3 - (int)no_space.size()), ' ');
  }

  // Compare two atoms for bond ordering: match MMDB behavior
  // Sort by element atomic number ascending, then by trimmed name length descending
  static bool bond_order_first(const gemmi::Atom *a1, const gemmi::Atom *a2) {
    int z1 = a1->element.atomic_number();
    int z2 = a2->element.atomic_number();
    if (z1 != z2) return z1 < z2;
    // Same element: longer name first (e.g., CA before C)
    std::string n1, n2;
    for (char c : a1->name) if (c != ' ') n1 += c;
    for (char c : a2->name) if (c != ' ') n2 += c;
    if (n1.size() != n2.size()) return n1.size() > n2.size();
    return n1 < n2;
  }

  void init_gemmi(const gemmi::Residue *residue_p) {
    filled_with_bond_order_data_only_flag = false;
    nuclear_distances_flag = false;

    if (!residue_p) return;

    // Collect non-hydrogen atoms
    std::vector<const gemmi::Atom*> atoms;
    for (const auto &atom : residue_p->atoms) {
      if (!atom.element.is_hydrogen()) {
        atoms.push_back(&atom);
      }
    }

    // Build bonds using proximity with known covalent distances
    std::vector<atom_pair_t> bond_pairs;
    for (size_t i = 0; i < atoms.size(); i++) {
      for (size_t j = i + 1; j < atoms.size(); j++) {
        const gemmi::Atom *a1 = atoms[i];
        const gemmi::Atom *a2 = atoms[j];

        clipper::Coord_orth p1(a1->pos.x, a1->pos.y, a1->pos.z);
        clipper::Coord_orth p2(a2->pos.x, a2->pos.y, a2->pos.z);
        double dist = std::sqrt((p1 - p2).lengthsq());

        // Use element-specific cutoffs for covalent bonds
        std::string e1 = a1->element.name();
        std::string e2 = a2->element.name();
        double cutoff = 2.0;

        if (dist < cutoff) {
          std::string type = "single";
          if (dist < 1.30 &&
              ((e1 == "C" && e2 == "O") || (e1 == "O" && e2 == "C")))
            type = "double";

          // Order atoms: first atom has lower atomic number, tie-break by name
          const gemmi::Atom *first = a1;
          const gemmi::Atom *second = a2;
          if (!bond_order_first(a1, a2)) { std::swap(first, second); }

          dict_bond_restraint_t br(pad_name(first->name), pad_name(second->name),
                                   type, dist, 0.02, 0, false, 0.0, false);
          bond_restraint.push_back(br);
          bond_pairs.push_back(atom_pair_t(first, second));
        }
      }
    }

    // Find angles: pairs of bonds that share an atom
    for (unsigned int ibp = 0; ibp < bond_pairs.size(); ibp++) {
      for (unsigned int jbp = ibp; jbp < bond_pairs.size(); jbp++) {
        if (ibp != jbp) {
          const gemmi::Atom *shared = bond_pairs[ibp].shared_atom(bond_pairs[jbp]);
          if (shared) {
            const gemmi::Atom *at_1 = bond_pairs[ibp].at_1;
            const gemmi::Atom *at_2 = bond_pairs[ibp].at_2;
            const gemmi::Atom *at_3 = bond_pairs[jbp].at_1;

            if (at_1 == shared) {
              at_1 = bond_pairs[ibp].at_2;
              at_2 = bond_pairs[ibp].at_1;
            }
            if (at_3 == shared)
              at_3 = bond_pairs[jbp].at_2;

            if (at_1 && at_2 && at_3) {
              clipper::Coord_orth p1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
              clipper::Coord_orth p2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
              clipper::Coord_orth p3(at_3->pos.x, at_3->pos.y, at_3->pos.z);

              double angle = clipper::Util::rad2d(clipper::Coord_orth::angle(p1, p2, p3));

              if (angle > 0.001) {
                dict_angle_restraint_t ar(pad_name(at_1->name), pad_name(at_2->name),
                                          pad_name(at_3->name), angle, 3.0);
                angle_restraint.push_back(ar);
              }
            }
          }
        }
      }
    }

    // Build residue_info
    std::string comp_id = residue_p->name;
    int n_all = static_cast<int>(residue_p->atoms.size());
    int n_non_H = 0;
    for (const auto &atom : residue_p->atoms) {
      if (!atom.element.is_hydrogen()) n_non_H++;
    }

    residue_info = dict_chem_comp_t(comp_id, comp_id, comp_id, "monomer", n_all, n_non_H, ".");

    // Fill atom_info (all atoms including H)
    for (const auto &atom : residue_p->atoms) {
      std::string ele = atom.element.name();
      clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
      dict_atom da(atom.name, atom.name, ele, "", std::make_pair(false, clipper::Coord_orth()));
      da.model_Cartn = std::make_pair(true, pos);
      atom_info.push_back(da);
    }
  }
};

} // namespace coot