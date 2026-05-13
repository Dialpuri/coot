#include "function.hh"
#include <stdexcept>
#include <gemmi/model.hpp>

namespace pli {

std::pair<clipper::Coord_orth, clipper::Coord_orth>
pi_stacking_container_t::get_ring_pi_centre_points_gemmi(
    const std::vector<std::string> &ring_atom_names,
    const gemmi::Residue &res_ref,
    const std::string &chain_name) const {

    // dummy points, overwritten.
    clipper::Coord_orth pt_1(0, 0, 0);
    clipper::Coord_orth pt_2(0, 0, 0);

    // fill this vector
    std::vector<clipper::Coord_orth> aromatic_plane_points;
    for (const auto& atom_name : ring_atom_names) {
        bool found = false;
        for (const auto& atom : res_ref.atoms) {
            if (atom.name == atom_name) {
                clipper::Coord_orth at_pt(atom.pos.x, atom.pos.y, atom.pos.z);
                aromatic_plane_points.push_back(at_pt);
                found = true;
                break;
            }
        }
        if (!found) {
            // Build residue spec message using gemmi data
            std::string mess = "Not all aromatic atoms were found in residue ";
            mess += chain_name;
            mess += " "; 
            mess += std::to_string(res_ref.seqid.num.value); 
            mess += " ";
            throw std::runtime_error(mess);
        }
    }

    // can throw an exception
    std::pair<clipper::Coord_orth, clipper::Coord_orth> p =
        ring_centre_and_normal(aromatic_plane_points);
    double plane_offset = 1; // Angstroms plane normal offset, Clark & Labute 2007
    pt_1 = p.first + plane_offset * p.second;
    pt_2 = p.first - plane_offset * p.second;
    if (0) { 
        std::cout << "   Centre: " << p.first.x() << " " << p.first.y() << " " << p.first.z() << std::endl;
        std::cout << "   Normal: " << p.second.x() << " " << p.second.y() << " " << p.second.z() << std::endl;
        std::cout << " Pi point: " << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << std::endl;
    }
    return std::pair<clipper::Coord_orth, clipper::Coord_orth>(pt_1, pt_2);
}

// Helper function - can throw an exception if not enough points found in pts.
std::pair<clipper::Coord_orth, clipper::Coord_orth>
pi_stacking_container_t::ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const {
    clipper::Coord_orth centre(0,0,0);
    clipper::Coord_orth normal(0,0,0);

    if (pts.size() < 3) {
        std::string mess = "in ring_centre_and_normal() not enough point to calculate lsq plane";
        throw std::runtime_error(mess);
    }

    double sum_x = 0, sum_y = 0, sum_z = 0;
    for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
        sum_x += pts[ipl].x();
        sum_y += pts[ipl].y();
        sum_z += pts[ipl].z();
    }

    double divisor = 1.0/double(pts.size());
    centre = clipper::Coord_orth(sum_x*divisor, sum_y*divisor, sum_z*divisor);
    clipper::Matrix<double> mat(3,3);
    for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
        mat(0,0) += (pts[ipl].x()-centre.x()) * (pts[ipl].x()-centre.x());
        mat(1,1) += (pts[ipl].y()-centre.y()) * (pts[ipl].y()-centre.y());
        mat(2,2) += (pts[ipl].z()-centre.z()) * (pts[ipl].z()-centre.z());
        mat(0,1) += (pts[ipl].x()-centre.x()) * (pts[ipl].y()-centre.y());
        mat(0,2) += (pts[ipl].x()-centre.x()) * (pts[ipl].z()-centre.z());
        mat(1,2) += (pts[ipl].y()-centre.y()) * (pts[ipl].z()-centre.z());
    }
    mat(1,0) = mat(0,1);
    mat(2,0) = mat(0,2);
    mat(2,1) = mat(1,2);
    std::vector<double> eigens = mat.eigen(false);

    int eigen_index = 0;
    if (eigens[1] < eigens[eigen_index])
        eigen_index = 1;
    if (eigens[2] < eigens[eigen_index])
        eigen_index = 2;

    clipper::Coord_orth eigen_vec(mat(0, eigen_index),
                                  mat(1, eigen_index),
                                  mat(2, eigen_index));

    double sum_sq = 1e-20;
    int n = 0;
    for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
        clipper::Coord_orth d = pts[ipl] - centre;
        double dd = d * eigen_vec;
        sum_sq += dd*dd;
        n++;
    }
    double sigma = sqrt(sum_sq/n);
    if (0) {
        std::cout << "   Plane fit sigma: " << sigma << std::endl;
    }
    if (sigma > 0.5) {
        std::string mess = "pi_stacking: plane fit not good";
        throw std::runtime_error(mess);
    }
    normal = eigen_vec;

    return std::make_pair(centre, normal);
}

} // namespace pli