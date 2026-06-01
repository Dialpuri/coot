#include "function.hh"

namespace coot {
namespace hole {

vdw_radii_result assign_vdw_radii_gemmi(
    gemmi::Model& model,
    const coot::protein_geometry& geom)
{
    int imol = 0;
    bool use_vdwH_flag = false;

    std::map<std::pair<std::string, std::string>, double> cached_radii;
    vdw_radii_result result;

    result.radii.reserve(model.chains.size());

    for (gemmi::Chain& chain : model.chains) {
        std::vector<std::vector<double>> chain_radii;
        chain_radii.reserve(chain.residues.size());

        for (gemmi::Residue& residue : chain.residues) {
            std::vector<double> res_radii;
            res_radii.reserve(residue.atoms.size());

            std::string residue_name = residue.name;

            for (gemmi::Atom& atom : residue.atoms) {
                std::string padded_atom_name = pad_atom_name(atom.name);

                // try cache first
                std::pair<std::string, std::string> p(padded_atom_name, residue_name);
                auto cit = cached_radii.find(p);

                double radius;
                if (cit != cached_radii.end()) {
                    radius = cit->second;
                } else {
                    radius = geom.get_vdw_radius(padded_atom_name, residue_name, imol, use_vdwH_flag);
                    cached_radii[p] = radius;
                }

                if (radius > 0) {
                    res_radii.push_back(radius);
                } else {
                    // make a reasonable default based on element
                    double default_radius = 1.7;
                    std::string ele = atom.element.name();
                    if (ele == "N")
                        default_radius = 1.55;
                    if (ele == "O")
                        default_radius = 1.52;
                    if (ele == "H")
                        default_radius = 1.2;
                    res_radii.push_back(default_radius);
                }
            }
            chain_radii.push_back(std::move(res_radii));
        }
        result.radii.push_back(std::move(chain_radii));
    }

    return result;
}

} // namespace hole
} // namespace coot