#include "function.hh"
#include <gemmi/pdb.hpp>
#include <iostream>

namespace coot {
namespace molecule_t {
namespace modification_info_t {

save_info_t::save_info_t(const std::string &file_name, const std::string &modification_info_string)
    : file_name(file_name), modification_info_string(modification_info_string) {
  structure = std::make_unique<gemmi::Structure>();
  try {
    *structure = gemmi::read_pdb_file(file_name);
  } catch (...) {
    structure.reset(); // leave as empty structure
  }
}

const gemmi::Structure* save_info_t::get_mol_gemmi() const {
  return structure ? structure.get() : nullptr;
}

} // namespace modification_info_t
} // namespace molecule_t
} // namespace coot