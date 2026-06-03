#include "function.hh"
#include <gemmi/math.hpp>
#include <gemmi/unitcell.hpp>

namespace coot {

std::vector<gemmi::Atom> trans_sel_gemmi(
    gemmi::Structure* st,
    std::vector<gemmi::Atom> atoms,
    int symm_no,
    int symm_x, int symm_y, int symm_z,
    int us, int vs, int ws)
{
    if (!st || st->models.empty()) {
        return atoms;
    }

    const gemmi::UnitCell& cell = st->cell;

    // Get symmetry operators from space group
    const gemmi::SpaceGroup* sg = st->find_spacegroup();
    if (!sg) {
        return atoms;
    }

    const gemmi::GroupOps& group_ops = sg->operations();
    const std::vector<gemmi::Op>& ops = group_ops.sym_ops;

    if (symm_no < 0 || symm_no >= static_cast<int>(ops.size())) {
        return atoms;
    }

    // Original code:
    //   GetTMatrix(my_matt, isym, x, y, z)  -- symmetry + lattice translation
    //   GetTMatrix(mol_to_origin, 0, -us, -vs, -ws)  -- lattice translation
    //   atom->Transform(mol_to_origin);
    //   atom->Transform(my_matt);
    //
    // This means: mol_to_origin first, then my_matt.
    // In fractional terms: frac -> frac + (-us,-vs,-ws) -> symop -> frac + (x,y,z) -> orthogonal

    gemmi::Op sym_op = ops[symm_no];

    // Create lattice translation for the symmetry (symm_x, symm_y, symm_z)
    gemmi::Op sym_latt;
    constexpr int DEN = gemmi::Op::DEN;
    sym_latt.rot = {{{DEN, 0, 0}, {0, DEN, 0}, {0, 0, DEN}}};
    sym_latt.tran = {symm_x * DEN, symm_y * DEN, symm_z * DEN};
    sym_latt.notation = ' ';

    // Combine: sym_op then lattice translation
    gemmi::Op sym_with_latt = sym_op.combine(sym_latt);

    // Create lattice translation for mol_to_origin (-us, -vs, -ws)
    gemmi::Op origin_trans;
    origin_trans.rot = {{{DEN, 0, 0}, {0, DEN, 0}, {0, 0, DEN}}};
    origin_trans.tran = {-us * DEN, -vs * DEN, -ws * DEN};
    origin_trans.notation = ' ';

    // Combined: origin_trans first, then sym_with_latt
    // apply origin_trans then sym_with_latt: origin_trans.combine(sym_with_latt)
    gemmi::Op combined = origin_trans.combine(sym_with_latt);

    // Apply combined transform to each atom
    for (gemmi::Atom& atom : atoms) {
        gemmi::Fractional frac = cell.fractionalize(atom.pos);
        std::array<double, 3> frac_in = {frac.x, frac.y, frac.z};
        std::array<double, 3> frac_sym = combined.apply_to_xyz(frac_in);
        atom.pos = cell.orthogonalize(
            gemmi::Fractional(frac_sym[0], frac_sym[1], frac_sym[2]));
    }

    return atoms;
}

} // namespace coot