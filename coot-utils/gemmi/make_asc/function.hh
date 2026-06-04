#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <unordered_map>

//! gemmi-compatible atom selection container (mirrors atom_selection_container_t)
struct atom_selection_container_t {
  const gemmi::Structure* mol;
  int n_selected_atoms;
  std::vector<int> atom_selection;
  std::string read_error_message;
  int read_success;
  int SelectionHandle;
  int UDDAtomIndexHandle;
  int UDDOldAtomIndexHandle;
};

atom_selection_container_t make_asc_gemmi(const gemmi::Structure& st, bool transfer_atom_index_flag);