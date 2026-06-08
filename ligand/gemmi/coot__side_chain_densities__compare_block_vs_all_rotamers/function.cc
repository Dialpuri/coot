#include "function.hh"

namespace coot {
namespace side_chain_densities {

std::map<std::string, std::pair<std::string, double>>
compare_block_vs_all_rotamers_gemmi(
    density_box_t block,
    gemmi::Residue *residue_p,
    const std::string &data_dir,
    const std::pair<bool, std::vector<std::pair<std::string, std::string>>> &rotamer_limits,
    const clipper::Xmap<float> &xmap) {

  std::map<std::string, std::pair<std::string, double>> probability_map;

  std::string glob_pattern = "*";
  std::vector<std::string> dirs = coot::util::glob_files(data_dir, glob_pattern);

  for (std::size_t idir = 0; idir < dirs.size(); idir++) {
    const std::string &res_dir = dirs[idir];

    std::vector<std::string> rot_dirs = coot::util::glob_files(res_dir, glob_pattern);
    for (std::size_t jdir = 0; jdir < rot_dirs.size(); jdir++) {
      const std::string &rot_dir = rot_dirs[jdir];

      std::string res = coot::util::file_name_non_directory(res_dir);
      std::string rot = coot::util::file_name_non_directory(rot_dir);
      std::string key = res + ":" + rot;

      bool do_it = false;

      if (!rotamer_limits.first) {
        // so, don't apply the limits
        do_it = true;
      } else {
        if (rotamer_limits.first) {
          // testing path
          if (!rotamer_limits.second.empty()) {
            bool found = false;
            for (std::size_t i = 0; i < rotamer_limits.second.size(); i++) {
              std::string limit_key =
                  rotamer_limits.second[i].first + ":" + rotamer_limits.second[i].second;
              if (limit_key == key) {
                found = true;
                break;
              }
            }
            if (found)
              do_it = true;
          }
        }
      }

      if (do_it)
        if (rot == "none")
          if (res != "GLY")
            if (res != "ALA")
              do_it = false;

      if (do_it) {

        std::pair<bool, double> p =
            coot::compare_block_vs_rotamer_gemmi(block, residue_p, rot_dir, xmap);
        if (p.first) {
          probability_map[key] = std::make_pair(res, p.second);
        }
      }
    }
  }
  return probability_map;
}

} // namespace side_chain_densities
} // namespace coot
