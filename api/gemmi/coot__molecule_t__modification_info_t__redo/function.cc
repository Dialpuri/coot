#include "function.hh"
#include <iostream>

namespace coot {
namespace molecule_t {
namespace modification_info_t {

const gemmi::Structure* redo_gemmi(
    std::vector<save_info_t>& save_info,
    int& modification_index) {
  
  std::cout << "coot::molecule_t::modification_info_t::redo()" << std::endl;
  
  int idx = modification_index + 1;
  if (idx > int(save_info.size()))
    idx = save_info.size();
  std::cout << ":::::::::::: in redo() modification_index: " << modification_index
            << " idx of molecule to change to: " << idx << std::endl;
  
  const gemmi::Structure* MMDBManager = nullptr;
  if (idx >= 0) {
    if (idx < int(save_info.size())) {
      MMDBManager = save_info[idx].get_mol_gemmi();
      modification_index = idx;
    }
  }
  return MMDBManager;
}

} // namespace modification_info_t
} // namespace molecule_t
} // namespace coot