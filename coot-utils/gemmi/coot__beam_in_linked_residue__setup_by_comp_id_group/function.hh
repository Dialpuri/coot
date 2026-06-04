#pragma once

#include <string>
#include <iostream>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_nth_residue/gemmi/function.hh"
#include <coot/utils/coot-utils.hh>

namespace coot {

namespace beam_in_linked_residue {

inline bool setup_by_comp_id_group_gemmi(
    const std::string& comp_id_ref,
    const std::string& group_new,
    const std::string& link_type,
    gemmi::CRA& template_res_ref,
    gemmi::CRA& template_res_mov,
    bool& have_template) {

  bool status = false;

  std::string file_name = comp_id_ref;
  file_name += "-";
  file_name += group_new;
  file_name += "-via-";
  file_name += link_type;
  file_name += ".pdb";

  std::string pkgdatadir = coot::package_data_dir();
  std::string full_path_pdb_filename = pkgdatadir + "/" + file_name;

  if (!coot::file_exists(full_path_pdb_filename)) {
    std::cout << "WARNING:: link template file " << full_path_pdb_filename
              << " does not exist " << std::endl;
  } else {
    gemmi::Structure t_st = gemmi::read_pdb_file(full_path_pdb_filename);

    template_res_ref = coot::util::get_nth_residue_gemmi(1, t_st);
    if (template_res_ref.atom == nullptr) {
      std::cout << "ERROR:: failed to find residue with comp_id " << comp_id_ref
                << " in " << full_path_pdb_filename << std::endl;
    } else {
      template_res_mov = coot::util::get_nth_residue_gemmi(2, t_st);
      if (template_res_mov.atom == nullptr) {
        std::cout << "ERROR:: failed to find (adding) residue with comp_id " << group_new
                  << " in " << full_path_pdb_filename << std::endl;
      } else {
        have_template = true;
        status = true;
      }
    }
  }

  return status;
}

} // namespace beam_in_linked_residue

} // namespace coot