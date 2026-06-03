#pragma once

#include "clipper/core/coords.h"
#include <gemmi/model.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>

clipper::Coord_orth
closest_approach_gemmi(const clipper::Coord_orth &moving_point,
                       const clipper::Coord_orth &reference_point,
                       const gemmi::Structure &mol) {
    clipper::Coord_orth pos(-1.0, -1.0, -1.0);
    double closest_dist_sq = 99999999999.9;

    const gemmi::UnitCell &cell = mol.cell;

    const gemmi::SpaceGroup *sg = mol.find_spacegroup();
    if (!sg) {
        return pos;
    }

    const gemmi::GroupOps &group_ops = sg->operations();
    const std::vector<gemmi::Op> &ops = group_ops.sym_ops;

    gemmi::Position moving_orth(moving_point.x(), moving_point.y(), moving_point.z());
    gemmi::Fractional moving_frac = cell.fractionalize(moving_orth);
    std::array<double, 3> frac_in = {moving_frac.x, moving_frac.y, moving_frac.z};

    for (int x_shift = -2; x_shift <= 2; x_shift++) {
        for (int y_shift = -2; y_shift <= 2; y_shift++) {
            for (int z_shift = -2; z_shift <= 2; z_shift++) {
                for (size_t ii = 0; ii < ops.size(); ii++) {
                    gemmi::Op sym_op = ops[ii];

                    // Create lattice translation
                    gemmi::Op trans;
                    constexpr int DEN = gemmi::Op::DEN;
                    trans.rot = {{{DEN, 0, 0}, {0, DEN, 0}, {0, 0, DEN}}};
                    trans.tran = {x_shift * DEN, y_shift * DEN, z_shift * DEN};
                    trans.notation = ' ';

                    gemmi::Op combined = trans.combine(sym_op);

                    std::array<double, 3> frac_sym = combined.apply_to_xyz(frac_in);
                    gemmi::Position trans_ort = cell.orthogonalize(
                        gemmi::Fractional(frac_sym[0], frac_sym[1], frac_sym[2]));

                    double dx = trans_ort.x - reference_point.x();
                    double dy = trans_ort.y - reference_point.y();
                    double dz = trans_ort.z - reference_point.z();
                    double d_sq = dx * dx + dy * dy + dz * dz;

                    if (d_sq < closest_dist_sq) {
                        closest_dist_sq = d_sq;
                        pos = clipper::Coord_orth(trans_ort.x, trans_ort.y, trans_ort.z);
                    }
                }
            }
        }
    }

    return pos;
}