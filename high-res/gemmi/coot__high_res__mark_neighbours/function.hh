#pragma once

#include <string>
#include <vector>

namespace coot {
namespace high_res {

/*  Recursively mark atoms (by flat index) whose atom name matches.
 *
 *  This is the gemmi-port of coot::high_res::mark_neighbours.
 *  Instead of MMDB UDData we use a plain std::vector<int> (group)
 *  indexed by flat-atom-position, initialized to -1.
 *
 *  Parameters
 *    iatom      – flat index into atom_names / group
 *    igroup     – group id to write into group[iatom]
 *    atom_name  – atom name to match (PDB-style padded, e.g. " N  ")
 *    neighbours – adjacency list  neighbours[i] = list of neighbor indices
 *    atom_names – flat list of atom names (one per atom, PDB-padded)
 *    group      – mutable group array, initialized to -1 before first call
 */
void mark_neighbours_gemmi(
    int iatom,
    int igroup,
    const std::string &atom_name,
    const std::vector<std::vector<int>> &neighbours,
    const std::vector<std::string> &atom_names,
    std::vector<int> &group);

} // namespace high_res
} // namespace coot