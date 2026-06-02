#include "function.hh"
#include <gemmi/math.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residue_spec_t__select_atoms/gemmi/function.hh"

namespace coot {

std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>
symmetry_contacts_gemmi(
    const residue_spec_t& res_central,
    gemmi::Structure* st,
    float d
) {
    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> v;

    if (!st || st->models.empty()) {
        return v;
    }

    gemmi::Model& model = st->models[0];
    const gemmi::UnitCell& cell = st->cell;

    // For non-crystal structures, MMDB's GetTMatrix fails, so return empty
    if (!cell.is_crystal()) {
        return v;
    }

    // Get symmetry operators from space group
    const gemmi::SpaceGroup* sg = st->find_spacegroup();
    if (!sg) {
        return v;
    }

    const gemmi::GroupOps& group_ops = sg->operations();
    const std::vector<gemmi::Op>& ops = group_ops.sym_ops;
    int n_symm = static_cast<int>(ops.size());
    if (n_symm == 0) {
        return v;
    }

    int shift_lim = 2;

    // Select atoms in the central residue using the verified gemmi port
    std::vector<gemmi::Atom*> res_atoms = coot::select_atoms_gemmi(res_central, st);
    if (res_atoms.empty()) {
        return v;
    }

    // Collect all atoms in the model
    std::vector<gemmi::Atom*> all_atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                all_atoms.push_back(&atom);
            }
        }
    }

    if (all_atoms.empty()) {
        return v;
    }

    // For each symmetry operation and lattice translation
    for (int x_shift = -shift_lim; x_shift <= shift_lim; x_shift++) {
        for (int y_shift = -shift_lim; y_shift <= shift_lim; y_shift++) {
            for (int z_shift = -shift_lim; z_shift <= shift_lim; z_shift++) {
                for (int i_symm = 0; i_symm < n_symm; i_symm++) {
                    // Skip identity (0,0,0,0)
                    if (!(x_shift == 0 && y_shift == 0 && z_shift == 0 && i_symm == 0)) {

                        gemmi::Op sym_op = ops[i_symm];

                        // Create a lattice translation op
                        gemmi::Op trans;
                        constexpr int DEN = gemmi::Op::DEN;
                        trans.rot = {{{DEN,0,0},{0,DEN,0},{0,0,DEN}}};
                        trans.tran = {x_shift * DEN, y_shift * DEN, z_shift * DEN};
                        trans.notation = ' ';

                        // Combine: translation ∘ sym_op (sym_op first, then translation)
                        gemmi::Op combined = trans.combine(sym_op);

                        // For each atom in the central residue, find contacts
                        // with all atoms transformed by the symmetry
                        for (gemmi::Atom* at1 : res_atoms) {
                            gemmi::Position pos1(at1->pos.x, at1->pos.y, at1->pos.z);

                            for (gemmi::Atom* at2 : all_atoms) {
                                // Convert to fractional, apply symmetry, convert back
                                gemmi::Fractional frac = cell.fractionalize(
                                    gemmi::Position(at2->pos.x, at2->pos.y, at2->pos.z));
                                std::array<double, 3> frac_in = {frac.x, frac.y, frac.z};
                                std::array<double, 3> frac_sym = combined.apply_to_xyz(frac_in);
                                gemmi::Position at2_sym = cell.orthogonalize(
                                    gemmi::Fractional(frac_sym[0], frac_sym[1], frac_sym[2]));

                                double dx = pos1.x - at2_sym.x;
                                double dy = pos1.y - at2_sym.y;
                                double dz = pos1.z - at2_sym.z;
                                double dist_sq = dx*dx + dy*dy + dz*dz;

                                if (dist_sq <= d * d && dist_sq > 0.0) {
                                    v.emplace_back(at1, at2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return v;
}

} // namespace coot