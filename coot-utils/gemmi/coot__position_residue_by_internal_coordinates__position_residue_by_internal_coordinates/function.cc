#include "function.hh"

#include <gemmi/math.hpp>
#include <cmath>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__position_residue_by_internal_coordinates__get_atom/gemmi/function.hh"

namespace coot {

position_residue_by_internal_coordinates_gemmi::position_residue_by_internal_coordinates_gemmi(
    gemmi::Residue *residue_ref,
    gemmi::Residue *residue_moving,
    const atom_name_quad &quad,
    const double &bond_length,
    const double &bond_angle,
    const double &bond_torsion)
    : residue_ref_(residue_ref)
    , residue_moving_(residue_moving)
    , quad_(quad)
    , bond_length_(bond_length)
    , bond_angle_(bond_angle)
    , bond_torsion_(bond_torsion)
    , at_0_(nullptr)
    , at_1_(nullptr)
    , at_2_(nullptr)
    , at_3_(nullptr)
{
    // first, can we find the atoms in the residues?
    if (residue_ref_) {
        if (residue_moving_) {
            at_0_ = coot::position_residue_by_internal_coordinates_get_atom_gemmi(
                residue_ref_, residue_moving_, quad_, 0);
            at_1_ = coot::position_residue_by_internal_coordinates_get_atom_gemmi(
                residue_ref_, residue_moving_, quad_, 1);
            at_2_ = coot::position_residue_by_internal_coordinates_get_atom_gemmi(
                residue_ref_, residue_moving_, quad_, 2);
            at_3_ = coot::position_residue_by_internal_coordinates_get_atom_gemmi(
                residue_ref_, residue_moving_, quad_, 3);

            if (at_0_ && at_1_ && at_2_ && at_3_) {
                // OK, we have found the atoms of the residue
            }
        }
    }
}

bool position_residue_by_internal_coordinates_gemmi::move_moving_residue() {
    // If atoms weren't found in constructor, return false
    if (!at_0_ || !at_1_ || !at_2_ || !at_3_) return false;
    if (!residue_ref_ || !residue_moving_) return false;

    // Convert degrees to radians
    double bond_angle_rad = bond_angle_ * M_PI / 180.0;
    double bond_torsion_rad = bond_torsion_ * M_PI / 180.0;

    // Reference frame from atoms 0,1,2
    gemmi::Position p0(at_0_->pos);
    gemmi::Position p1(at_1_->pos);
    gemmi::Position p2(at_2_->pos);

    // v1 = p2 - p1 (z-axis direction)
    gemmi::Vec3 v1 = p2 - p1;
    double d1 = v1.length();
    if (d1 < 1e-8) return false;

    // v2 = p0 - p1
    gemmi::Vec3 v2 = p0 - p1;
    gemmi::Vec3 n1 = v1 / d1;  // z-axis

    // x-axis: component of v2 perpendicular to v1
    gemmi::Vec3 n2 = v2 - v1.dot(n1) * n1;
    double d2 = n2.length();
    if (d2 < 1e-8) return false;
    n2 = n2 / d2;

    // y-axis: cross product
    gemmi::Vec3 n3 = n1.cross(n2);

    // Compute target position for atom 3
    double r  = bond_length_;
    double theta = bond_angle_rad;
    double phi = bond_torsion_rad;

    double dx =  r * sin(theta) * cos(phi);
    double dy =  r * sin(theta) * sin(phi);
    double dz = -r * cos(theta);

    gemmi::Vec3 disp(dx * n2.x + dy * n3.x + dz * n1.x,
                     dx * n2.y + dy * n3.y + dz * n1.y,
                     dx * n2.z + dy * n3.z + dz * n1.z);

    gemmi::Position new_p3(p2.x + disp.x, p2.y + disp.y, p2.z + disp.z);

    // Compute rotation that maps old local frame to new
    gemmi::Vec3 new_v2 = new_p3 - p2;
    double new_d = new_v2.length();
    if (new_d < 1e-8) return false;

    gemmi::Vec3 new_n2 = new_v2 / new_d;
    gemmi::Vec3 new_n3 = n1.cross(new_n2);
    double new_d3 = new_n3.length();
    if (new_d3 < 1e-8) return false;
    new_n3 = new_n3 / new_d3;

    // Rotate all atoms in the moving residue around p2
    for (gemmi::Atom& a : residue_moving_->atoms) {
        gemmi::Vec3 rel = a.pos - p2;
        // Project into old local frame
        double cx = n2.dot(rel);
        double cy = n3.dot(rel);
        double cz = n1.dot(rel);
        // Express in new local frame
        gemmi::Vec3 rotated(cx * new_n2.x + cy * new_n3.x + cz * n1.x,
                            cx * new_n2.y + cy * new_n3.y + cz * n1.y,
                            cx * new_n2.z + cy * new_n3.z + cz * n1.z);
        a.pos = gemmi::Position(p2.x + rotated.x, p2.y + rotated.y, p2.z + rotated.z);
    }

    return true;
}

} // namespace coot