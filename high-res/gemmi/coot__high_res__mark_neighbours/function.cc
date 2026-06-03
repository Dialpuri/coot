#include "function.hh"

namespace coot {
namespace high_res {

void mark_neighbours_gemmi(
    int iatom,
    int igroup,
    const std::string &atom_name,
    const std::vector<std::vector<int>> &neighbours,
    const std::vector<std::string> &atom_names,
    std::vector<int> &group) {

  int ig = group[iatom];
  if (ig == -1) {
    if (atom_names[iatom] == atom_name) {
      group[iatom] = igroup;
      const std::vector<int> &n = neighbours[iatom];
      for (int in = 0; in < static_cast<int>(n.size()); in++) {
        mark_neighbours_gemmi(n[in], igroup, atom_name, neighbours, atom_names, group);
      }
    }
  }
}

} // namespace high_res
} // namespace coot