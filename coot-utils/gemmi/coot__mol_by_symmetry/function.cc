#include "function.hh"

namespace coot {

gemmi::Structure mol_by_symmetry_gemmi(
    const gemmi::Structure& mol,
    gemmi::UnitCell cell,
    gemmi::FTransform ftrans,
    std::vector<int> pre_shift_to_origin_abc)
{
    // Deep copy the structure
    gemmi::Structure mol2(mol);

    // Compute origin shift from pre_shift_to_origin_abc (lattice translation)
    gemmi::Vec3 origin_shift_orth(0, 0, 0);
    if (pre_shift_to_origin_abc.size() == 3) {
        gemmi::Position p = cell.orthogonalize(
            gemmi::Fractional(static_cast<double>(pre_shift_to_origin_abc[0]),
                              static_cast<double>(pre_shift_to_origin_abc[1]),
                              static_cast<double>(pre_shift_to_origin_abc[2])));
        origin_shift_orth = gemmi::Vec3(p.x, p.y, p.z);
    }

    // Convert fractional transform to orthogonal transform
    gemmi::Transform rtop = cell.orthogonalize_transform(ftrans);

    // Apply transformation to every atom in every model
    for (gemmi::Model& model : mol2.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    gemmi::Vec3 co(atom.pos);
                    co -= origin_shift_orth;           // shift to origin
                    co = rtop.apply(co);                    // apply symmetry transform
                    co += origin_shift_orth;           // shift back
                    atom.pos = gemmi::Position(co);
                }
            }
        }
    }

    return mol2;
}

} // namespace coot