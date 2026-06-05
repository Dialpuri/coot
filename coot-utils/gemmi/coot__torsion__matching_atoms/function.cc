#include "function.hh"

namespace coot {
namespace torsion {

std::vector<gemmi::Atom*> matching_atoms_gemmi(
    gemmi::Chain* chain,
    gemmi::Residue* residue,
    const std::pair<int, coot::atom_spec_t>& atom_1,
    const std::pair<int, coot::atom_spec_t>& atom_2,
    const std::pair<int, coot::atom_spec_t>& atom_3,
    const std::pair<int, coot::atom_spec_t>& atom_4) {

    std::vector<gemmi::Atom*> v;

    gemmi::Atom* catom_1 = nullptr;
    gemmi::Atom* catom_2 = nullptr;
    gemmi::Atom* catom_3 = nullptr;
    gemmi::Atom* catom_4 = nullptr;

    for (gemmi::Atom& atom : residue->atoms) {
        gemmi::CRA cra{chain, residue, &atom};
        if (atom_1.second.matches_spec_gemmi(cra))
            catom_1 = &atom;
        if (atom_2.second.matches_spec_gemmi(cra))
            catom_2 = &atom;
        if (atom_3.second.matches_spec_gemmi(cra))
            catom_3 = &atom;
        if (atom_4.second.matches_spec_gemmi(cra))
            catom_4 = &atom;
    }

    if (catom_1 && catom_2 && catom_3 && catom_4) {
        v.push_back(catom_1);
        v.push_back(catom_2);
        v.push_back(catom_3);
        v.push_back(catom_4);
    }

    return v;
}

} // namespace torsion
} // namespace coot