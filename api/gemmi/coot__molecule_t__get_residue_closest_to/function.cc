#include "function.hh"
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

double distance_sq(const clipper::Coord_orth& co, const gemmi::Atom& atom) {
    double dx = atom.pos.x - co.x();
    double dy = atom.pos.y - co.y();
    double dz = atom.pos.z - co.z();
    return dx * dx + dy * dy + dz * dz;
}

residue_spec_t molecule_t::get_residue_closest_to_gemmi(const gemmi::Structure& st,
                                                        const clipper::Coord_orth& co) const {
    residue_spec_t spec;

    double d_best = 999999999999.0;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                for (const gemmi::Atom& atom : residue.atoms) {
                    // Skip terminal atoms (Ter) - in MMDB, isTer() checks for termination
                    // In gemmi, we check if atom name is "TER"
                    bool is_ter = (atom.name == "TER");
                    if (!is_ter) {
                        double dd = distance_sq(co, atom);
                        if (dd < d_best) {
                            d_best = dd;
                            spec = residue_spec_t(&chain, &residue);
                        }
                    }
                }
            }
        }
    }

    return spec;
}

} // namespace coot