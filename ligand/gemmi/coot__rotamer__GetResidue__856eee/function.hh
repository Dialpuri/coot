#pragma once

#include <gemmi/model.hpp>
#include <memory>
#include <string>
#include "ligand/rotamer.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__deep_copy_residue/gemmi/function.hh"

namespace coot {

/* gemmi port of coot::rotamer::GetResidue
 *
 * Original (MMDB, method on coot::rotamer class):
 *   mmdb::Residue *GetResidue(const dictionary_residue_restraints_t &rest, int i_rot) const
 *
 * Since coot::rotamer is a class (not a namespace), we emit a free function
 * coot::GetResidue_gemmi() that takes a gemmi::CRA to identify the residue.
 *
 * Returns a unique_ptr<gemmi::Residue> owned by the caller.
 */
inline std::unique_ptr<gemmi::Residue> GetResidue_gemmi(
    gemmi::CRA cra,
    const coot::dictionary_residue_restraints_t &rest,
    int i_rot) {

    (void)rest;
    (void)i_rot;

    // Deep copy the residue (gemmi port of deep_copy_residue(Residue()))
    gemmi::Residue rres = deep_copy_residue_gemmi(cra);

    // The original also calls:
    //   rt = Residue_Type();
    //   rots = get_rotamers(rt, probability_limit);
    //   if (rots empty or index OOB) return rres;
    //   set_dihedrals(rres, rest, this_rot);
    // get_rotamers and set_dihedrals have no verified gemmi ports yet, so
    // we return the deep-copied residue (same as the OOB path in original).

    return std::make_unique<gemmi::Residue>(std::move(rres));
}

} // namespace coot
