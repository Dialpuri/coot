#include "function.hh"

namespace coot {

std::vector<gemmi::Residue*>
residues_near_position(const clipper::Coord_orth &pt,
                       const gemmi::Structure &st,
                       double radius) {
    std::vector<gemmi::Residue*> v;
    
    if (st.models.empty())
        return v;
    
    gemmi::Model& model = const_cast<gemmi::Model&>(st.models[0]);
    
    for (auto& chain : model.chains) {
        for (auto& residue : chain.residues) {
            for (auto& atom : residue.atoms) {
                clipper::Coord_orth at_pt(atom.pos.x, atom.pos.y, atom.pos.z);
                double d = clipper::Coord_orth::length(pt, at_pt);
                if (d < radius) {
                    v.push_back(&residue);
                    break;
                }
            }
        }
    }
    
    return v;
}

} // namespace coot