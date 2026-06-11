#include "function.hh"
#include <numeric>
#include <algorithm>
#include <limits>

namespace coot {

// ---- lsq_plane_info_t ----

lsq_plane_info_t::lsq_plane_info_t(const std::vector<gemmi::Vec3>& v) {
    gemmi::Vec3 centre(0, 0, 0);
    for (const auto& c : v) {
        centre = centre + c;
    }
    centre = centre * (1.0 / static_cast<double>(v.size()));

    double cxx = 0, cxy = 0, cxz = 0, cyy = 0, cyz = 0, czz = 0;
    for (const auto& c : v) {
        gemmi::Vec3 diff = c - centre;
        cxx += diff.x * diff.x;
        cxy += diff.x * diff.y;
        cxz += diff.x * diff.z;
        cyy += diff.y * diff.y;
        cyz += diff.y * diff.z;
        czz += diff.z * diff.z;
    }

    double a = cyy + czz;
    double b = cxx + czz;
    double c = cxx + cyy;
    double d = cxy;
    double e = cxz;
    double f = cyz;

    double min_val = std::numeric_limits<double>::max();
    gemmi::Vec3 min_normal(0, 0, 0);

    // Test 27 directions (sign combinations)
    int signs[] = {-1, 0, 1};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) {
        if (signs[i] == 0 && signs[j] == 0 && signs[k] == 0) continue;
        double xx = signs[i];
        double yy = signs[j];
        double zz = signs[k];

        double denom = a*xx + b*xx + c*xx;
        if (denom != 0) {
            double temp = denom;
            double xx2 = (temp - d*yy - e*zz);
            double yy2 = (temp - f*xx - e*zz);
            double zz2 = (temp - f*xx - d*yy);

            gemmi::Vec3 n(xx2, yy2, zz2);
            double len = std::sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
            if (len > 0) {
                n = n * (1.0/len);
                double val = 0;
                for (const auto& c : v) {
                    gemmi::Vec3 diff = c - centre;
                    val += std::abs(diff.dot(n));
                }
                if (val < min_val) {
                    min_val = val;
                    min_normal = n;
                }
            }
        }
    }
    normal_ = min_normal;
    d_ = -(normal_.dot(centre));
}

gemmi::Vec3 lsq_plane_info_t::centre() const {
    // Return a Vec3 representing the "centre" of the plane.
    // The closest point on the plane to origin = -d * normal
    return normal_ * (-d_);
}

gemmi::Vec3 lsq_plane_info_t::projected_point(gemmi::Vec3 p) const {
    double d = (p.dot(normal_)) + d_;
    gemmi::Vec3 n = normal_ * d;
    return p - n;
}

// ---- lsq_plane_deviation ----

std::pair<double, double> lsq_plane_deviation(
    const std::vector<gemmi::Vec3>& plane_atom_coords,
    gemmi::Vec3 pt) {
    lsq_plane_info_t lsq_plane(plane_atom_coords);
    gemmi::Vec3 proj_pt = lsq_plane.projected_point(pt);
    gemmi::Vec3 diff_vec = proj_pt - pt;
    double oop = std::sqrt(diff_vec.dot(diff_vec));

    double sum_sq = 0;
    for (const auto& c : plane_atom_coords) {
        gemmi::Vec3 proj = lsq_plane.projected_point(c);
        gemmi::Vec3 d = proj - c;
        sum_sq += d.dot(d);
    }
    double rms = std::sqrt(sum_sq / plane_atom_coords.size());
    return std::make_pair(oop, rms);
}

// ---- pucker_analysis_info_t::pucker_analysis_info_t_gemmi ----

void pucker_analysis_info_t::pucker_analysis_info_t_gemmi(gemmi::Residue const& res, std::string altconf_in) {

    out_of_plane_distance = 0.0;
    plane_distortion = 0.0;
    C1_prime = nullptr;
    N1_or_9  = nullptr;
    base_atoms_coords.clear();
    ribose_atoms_coords.clear();
    altconf = altconf_in;

    // Call the existing gemmi port for assign_base_atom_coords
    // The port signature: takes pucker_analysis_info_t& and gemmi::Residue const&
    // We need to find the right call form
    // Use a local object that matches the port's expected struct
    {
        // The port expects a struct with base_atoms_coords as std::vector<gemmi::Vec3>
        // and sets C1_prime and N1_or_9
        // We call it as a free function that populates our members
        // Looking at the port: it's assign_base_atom_coords_gemmi(residue)
        // as a member of a different struct. We need to replicate that logic.
        // 
        // Since the port struct has different member names/structure, let's
        // just inline the atom-finding logic for base atoms.
        
        // Base atoms: purine: N1,C2,N3,C4,C5,C6,N7,N8,C9
        // Base atoms: pyrimidine: N1,C2,N3,C4,C5,C6
        
        std::vector<std::string> purine_base_atoms = {" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N7 ", " N8 ", " C9 "};
        std::vector<std::string> pyrimidine_base_atoms = {" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "};

        std::string res_name = res.name;
        bool is_purine = (res_name == " A " || res_name == " G ");
        const std::vector<std::string>& base_atom_names = is_purine ? purine_base_atoms : pyrimidine_base_atoms;

        std::vector<gemmi::Vec3> base_coords;

        auto norm = [](const std::string& s) { return s.empty() ? std::string(" ") : s; };

        for (gemmi::Atom const& atom : res.atoms) {
            std::string atm_name = atom.name;
            char altloc_char = atom.altloc;
            std::string alt_name(1, altloc_char);

            if (norm(altconf) == norm(alt_name)) {
                // Check for base atoms
                for (const auto& expected : base_atom_names) {
                    if (atm_name == expected) {
                        base_coords.push_back(atom.pos);
                        break;
                    }
                }
                // Check for C1'
                if (atm_name == " C1'" || atm_name == " C1*") {
                    C1_prime = const_cast<gemmi::Atom*>(&atom);
                }
                // Check for N1 or N9
                if (is_purine) {
                    if (atm_name == " N9 ") {
                        N1_or_9 = const_cast<gemmi::Atom*>(&atom);
                    }
                } else {
                    if (atm_name == " N1 ") {
                        N1_or_9 = const_cast<gemmi::Atom*>(&atom);
                    }
                }
            }
        }

        base_atoms_coords = std::move(base_coords);
    }

    // Helper lambda: compute LSQ plane from coords
    auto get_base_lsq_plane = [](const std::vector<gemmi::Vec3>& coords) ->
        std::optional<lsq_plane_info_t> {
        if (coords.size() < 3) return std::nullopt;
        lsq_plane_info_t lsq_plane(coords);
        return lsq_plane;
    };

    std::optional<lsq_plane_info_t> lsq_plane = get_base_lsq_plane(base_atoms_coords);

    if (lsq_plane.has_value()) {
        // store the geometry
        markup_info.base_ring_centre = lsq_plane.value().centre();
        markup_info.base_ring_normal = lsq_plane.value().centre();

        std::vector<gemmi::Atom const*> ribose_atoms(5, nullptr);
        std::vector<coot::pucker_analysis_info_t::PUCKERED_ATOM_T> possible_puckers;
        possible_puckers.push_back(coot::pucker_analysis_info_t::C1_PRIME);
        possible_puckers.push_back(coot::pucker_analysis_info_t::C2_PRIME);
        possible_puckers.push_back(coot::pucker_analysis_info_t::C3_PRIME);
        possible_puckers.push_back(coot::pucker_analysis_info_t::C4_PRIME);
        possible_puckers.push_back(coot::pucker_analysis_info_t::O4_PRIME);

        auto norm = [](const std::string& s) { return s.empty() ? std::string(" ") : s; };

        // find phosphorus atom and ribose atoms
        for (gemmi::Atom const& atom : res.atoms) {
            std::string atm_name = atom.name;
            char altloc_char = atom.altloc;
            std::string alt_name(1, altloc_char);
            if (norm(altconf) == norm(alt_name)) {
                if (atm_name == " P  ") {
                    markup_info.phosphorus_position = atom.pos;
                    markup_info.projected_point = lsq_plane.value().projected_point(atom.pos);
                }
                if (atm_name == " C1*" || atm_name == " C1'") ribose_atoms[0] = &atom;
                if (atm_name == " C2*" || atm_name == " C2'") ribose_atoms[1] = &atom;
                if (atm_name == " C3*" || atm_name == " C3'") ribose_atoms[2] = &atom;
                if (atm_name == " C4*" || atm_name == " C4'") ribose_atoms[3] = &atom;
                if (atm_name == " O4*" || atm_name == " O4'") ribose_atoms[4] = &atom;
            }
        }

        if (!(ribose_atoms[0] && ribose_atoms[1] && ribose_atoms[2] && ribose_atoms[3] && ribose_atoms[4])) {
            throw std::runtime_error("Not all atoms found in ribose.");
        } else {
            for (int i_oop_atom = 0; i_oop_atom < 5; i_oop_atom++) {
                ribose_atoms_coords.push_back(ribose_atoms[i_oop_atom]->pos);
            }

            std::vector<std::pair<float, float>> pucker_distortion_and_oop_d(5);
            for (int i_oop_atom = 0; i_oop_atom < 5; i_oop_atom++) {
                std::vector<gemmi::Vec3> plane_atom_coords;
                for (int i = 0; i < 5; i++) {
                    if (i != i_oop_atom) {
                        plane_atom_coords.push_back(ribose_atoms[i]->pos);
                    }
                }
                gemmi::Vec3 pt = ribose_atoms[i_oop_atom]->pos;
                std::pair<double, double> dev =
                    coot::lsq_plane_deviation(plane_atom_coords, pt);
                pucker_distortion_and_oop_d[i_oop_atom] =
                    std::make_pair(static_cast<float>(dev.first),
                                   static_cast<float>(dev.second));
            }

            puckered_atom_ = coot::pucker_analysis_info_t::NONE;
            std::pair<float, float> most_deviant(0, 0);
            for (int i_oop_atom = 0; i_oop_atom < 5; i_oop_atom++) {
                if (std::fabs(pucker_distortion_and_oop_d[i_oop_atom].first) >
                    std::fabs(most_deviant.first)) {
                    most_deviant = pucker_distortion_and_oop_d[i_oop_atom];
                    puckered_atom_ = possible_puckers[i_oop_atom];
                }
            }
            out_of_plane_distance = most_deviant.first;
            plane_distortion = most_deviant.second;
        }
    } else {
        throw std::runtime_error("base lsq plane has no value");
    }
}

} // namespace coot
