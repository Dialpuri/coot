#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__high_res__get_middle_pos/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__high_res__fill_globular_protein/gemmi/function.hh"

namespace coot {
namespace high_res {

coot::minimol::molecule high_res_gemmi(
    const coot::minimol::molecule &mol,
    const clipper::Coord_orth &given_centre,
    const gemmi::Structure &structure)
{
    coot::minimol::molecule globular_mol;

    if (mol.get_number_of_atoms() == 0) {
        return globular_mol;
    }

    auto middle_pos = coot::high_res::get_middle_pos_gemmi(structure);
    clipper::Coord_orth middle = middle_pos.first;

    coot_high_res_fill_globular_protein_gemmi(
        globular_mol, mol, given_centre, structure);

    return globular_mol;
}

} // namespace high_res
} // namespace coot
