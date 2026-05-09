#pragma once
#include <memory>
#include <string>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {
namespace modification_info_t {

class save_info_t {
public:
  save_info_t(const std::string &file_name, const std::string &modification_info_string);
  const gemmi::Structure* get_mol_gemmi() const;
  std::string file_name;
  std::string modification_info_string;
  std::unique_ptr<gemmi::Structure> structure;
};

} // namespace modification_info_t
} // namespace molecule_t
} // namespace coot