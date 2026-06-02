#pragma once
#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "utils/coot-utils.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_nth_residue/gemmi/function.hh"

namespace coot {
namespace beam_in_linked_residue {

inline bool setup_by_group_group_gemmi(
    const std::string& group_ref,
    const std::string& group_new,
    const std::string& link_type,
    gemmi::CRA& template_res_ref,
    gemmi::CRA& template_res_mov,
    bool& have_template)
{
    bool status = false;

    std::string file_name = group_ref;
    file_name += "-";
    file_name += group_new;
    file_name += "-via-";
    file_name += link_type;
    file_name += ".pdb";

    std::string pkgdatadir = coot::package_data_dir();
    std::string full_path_pdb_filename = pkgdatadir + "/" + file_name;

    std::cout << "debug:: setup_by_group() full_path_pdb_filename "
              << full_path_pdb_filename << std::endl;

    if (!coot::file_exists(full_path_pdb_filename)) {
        std::cout << "WARNING:: link template file " << full_path_pdb_filename
                  << " does not exist" << std::endl;
    } else {
        gemmi::Structure st = gemmi::read_pdb_file(full_path_pdb_filename);
        if (st.models.empty()) {
            std::cout << "ERROR:: on reading " << full_path_pdb_filename << std::endl;
        } else {
            template_res_ref = coot::util::get_nth_residue_gemmi(1, st);
            if (template_res_ref.residue == nullptr) {
                std::cout << "ERROR:: failed to find residue 1 in "
                          << full_path_pdb_filename << std::endl;
            } else {
                template_res_mov = coot::util::get_nth_residue_gemmi(2, st);
                if (template_res_mov.residue == nullptr) {
                    std::cout << "ERROR:: failed to find residue 2 in "
                              << full_path_pdb_filename << std::endl;
                } else {
                    have_template = true;
                    status = true;
                }
            }
        }
    }

    return status;
}

} // namespace beam_in_linked_residue
} // namespace coot