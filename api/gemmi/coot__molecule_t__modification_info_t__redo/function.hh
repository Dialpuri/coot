#pragma once
#include <string>
#include <vector>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <memory>

namespace coot {
namespace molecule_t {
namespace modification_info_t {

class save_info_t {
public:
  save_info_t(const std::string &file_name, const std::string &modification_info_string)
    : file_name(file_name), modification_info_string(modification_info_string) {}
  const gemmi::Structure* get_mol_gemmi() const {
    return structure.get();
  }
  std::string file_name;
  std::string modification_info_string;
  std::unique_ptr<gemmi::Structure> structure;
};

// Port of redo() to gemmi - returns the structure at modification_index+1
// Returns nullptr if no redo available (already at end of history)
const gemmi::Structure* redo_gemmi(
    std::vector<save_info_t>& save_info,
    int& modification_index);

} // namespace modification_info_t
} // namespace molecule_t
} // namespace coot