#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <clipper/core/clipper_util.h>
#include <string>
#include <vector>
#include <utility>
#include <map>

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
};

struct dict_angle_restraint_t {
  std::string atom_id_1_;
  std::string atom_id_2_;
  std::string atom_id_3_;
  double angle_;
  double angle_esd_;

  dict_angle_restraint_t() : angle_(0), angle_esd_(0) {}
};

class dictionary_residue_restraints_t {
public:
  std::vector<dict_bond_restraint_t> bond_restraint;
  std::vector<dict_angle_restraint_t> angle_restraint;
  std::vector<dict_atom> atom_info;
  dict_chem_comp_t residue_info;
  int has_partial_charges_flag;
  int read_number;
  bool filled_with_bond_order_data_only_flag;
  bool nuclear_distances_flag;

  dictionary_residue_restraints_t()
    : has_partial_charges_flag(0), read_number(0),
      filled_with_bond_order_data_only_flag(false), nuclear_distances_flag(false) {}

  std::string comp_id() const { return residue_info.comp_id; }
  int is_filled() const { return filled_with_bond_order_data_only_flag ? 1 : 0; }
  unsigned int number_of_atoms() const { return atom_info.size(); }

  // gemmi port of the (comp_id, read_number) constructor
  static dictionary_residue_restraints_t dictionary_residue_restraints_t_gemmi(
      const std::string &comp_id_in, int read_number_in) {
    dictionary_residue_restraints_t rest;
    rest.residue_info.comp_id = comp_id_in;
    rest.has_partial_charges_flag = 0;
    rest.read_number = read_number_in;
    rest.filled_with_bond_order_data_only_flag = false;
    rest.nuclear_distances_flag = false;
    return rest;
  }
};

} // namespace coot
