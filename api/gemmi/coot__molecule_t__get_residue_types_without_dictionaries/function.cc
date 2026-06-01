#include "function.hh"

namespace coot { namespace molecule_t {

std::vector<std::string> get_residue_types_without_dictionaries_gemmi(
    const gemmi::Model &model,
    const coot::protein_geometry &geom,
    int imol_no)
{
    std::vector<std::string> v;

    for (const gemmi::Chain &chain : model.chains) {
        for (const gemmi::Residue &residue : chain.residues) {
            std::string rn(residue.name);
            if (!geom.have_dictionary_for_residue_type_no_dynamic_add(rn, imol_no)) {
                v.push_back(rn);
            }
        }
    }
    return v;
}

}} // namespace coot::molecule_t