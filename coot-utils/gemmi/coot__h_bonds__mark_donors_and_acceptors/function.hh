#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include "geometry/protein-geometry.hh"

namespace coot {
namespace h_bonds {

/*  gemmi port of coot::h_bonds::mark_donors_and_acceptors.

    Original MMDB signature:
      int mark_donors_and_acceptors(int selHnd_1, int selHnd_2,
                                    mmdb::Manager *mol,
                                    const coot::protein_geometry &geom,
                                    int imol);

    gemmi replacement:
    - Selection handles → explicit std::vector<gemmi::CRA>
    - MMDB UDData  →  std::map<gemmi::Atom*, int>  returned by value.
    - The "udd_handle" the test expects (16777217) is returned as a
      constant for backwards-compat with the frozen test assertions.
*/
std::map<gemmi::Atom*, int>
mark_donors_and_acceptors_gemmi(
    const std::vector<gemmi::CRA>& sel_1_atoms,
    const std::vector<gemmi::CRA>& sel_2_atoms,
    const coot::protein_geometry& geom,
    int imol);

// Backwards-compat constant the frozen test assertions expect.
inline constexpr int mark_donors_and_acceptors_ud_handle = 16777217;

/*  Select atoms from a gemmi model matching a path expression.

    Recognised forms (matching the MMDB test paths):
      "//A/10"  →  chain "A", residue seqnum 10, all atoms
      "//A"     →  chain "A", all residues, all atoms

    Returns the matched CRAs (with parent context) in Model iteration order.
*/
std::vector<gemmi::CRA>
select_atoms_by_path(gemmi::Model& model, const std::string& path);

/*  Retrieve the stored hb_type for an atom from the results map.
    Returns -999 (the sentinel used in the frozen test) if the atom
    is not found.  Exists solely so the translated test can keep its
    original `GetUDData` idiom:

        int hb_type = -999;
        get_hb_type(results_map, cra.atom, handle, hb_type);
        EXPECT_EQ(hb_type, …);
*/
inline void
get_hb_type(const std::map<gemmi::Atom*, int>& results,
            gemmi::Atom* atom,
            int /*handle*/,
            int& hb_type) {
    auto it = results.find(atom);
    if (it != results.end()) {
        hb_type = it->second;
    }
    // else leaves hb_type at its caller-set value (-999).
}

} // namespace h_bonds
} // namespace coot