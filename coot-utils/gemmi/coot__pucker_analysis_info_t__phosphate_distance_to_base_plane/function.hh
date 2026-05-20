#pragma once
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

struct pucker_analysis_info_t {
    std::vector<gemmi::Vec3> base_atoms_coords;
    char altconf;

    pucker_analysis_info_t(gemmi::Residue* res, const std::string& alt)
        : altconf(alt.empty() ? ' ' : alt[0]) {
        assign_base_atom_coords_gemmi(res);
    }

    void assign_base_atom_coords_gemmi(gemmi::Residue* residue) {
        std::vector<std::string> cytidine_base_names;
        std::vector<std::string> uracil_base_names;
        std::vector<std::string> adenine_base_names;
        std::vector<std::string> guanine_base_names;

        cytidine_base_names.push_back(" N1 ");
        cytidine_base_names.push_back(" C2 ");
        cytidine_base_names.push_back(" N3 ");
        cytidine_base_names.push_back(" C4 ");
        cytidine_base_names.push_back(" C5 ");
        cytidine_base_names.push_back(" C6 ");
        cytidine_base_names.push_back(" O2 ");
        cytidine_base_names.push_back(" N4 ");

        uracil_base_names.push_back(" N1 ");
        uracil_base_names.push_back(" C2 ");
        uracil_base_names.push_back(" N3 ");
        uracil_base_names.push_back(" C4 ");
        uracil_base_names.push_back(" C5 ");
        uracil_base_names.push_back(" C6 ");
        uracil_base_names.push_back(" O2 ");
        uracil_base_names.push_back(" O4 ");

        adenine_base_names.push_back(" N9 ");
        adenine_base_names.push_back(" C8 ");
        adenine_base_names.push_back(" N7 ");
        adenine_base_names.push_back(" C5 ");
        adenine_base_names.push_back(" C4 ");
        adenine_base_names.push_back(" N1 ");
        adenine_base_names.push_back(" C2 ");
        adenine_base_names.push_back(" N3 ");
        adenine_base_names.push_back(" C6 ");
        adenine_base_names.push_back(" N6 ");

        guanine_base_names.push_back(" N9 ");
        guanine_base_names.push_back(" C8 ");
        guanine_base_names.push_back(" N7 ");
        guanine_base_names.push_back(" C5 ");
        guanine_base_names.push_back(" C4 ");
        guanine_base_names.push_back(" N1 ");
        guanine_base_names.push_back(" C2 ");
        guanine_base_names.push_back(" N3 ");
        guanine_base_names.push_back(" C6 ");
        guanine_base_names.push_back(" O6 ");
        guanine_base_names.push_back(" N2 ");

        base_atoms_coords.clear();

        std::vector<std::string> base_names;
        std::string residue_name = residue->name;

        if (residue_name == "C") base_names = cytidine_base_names;
        if (residue_name == "U") base_names = uracil_base_names;
        if (residue_name == "A") base_names = adenine_base_names;
        if (residue_name == "G") base_names = guanine_base_names;
        if (residue_name == "Cr") base_names = cytidine_base_names;
        if (residue_name == "Ur") base_names = uracil_base_names;
        if (residue_name == "Ar") base_names = adenine_base_names;
        if (residue_name == "Gr") base_names = guanine_base_names;
        if (residue_name == "CYT") base_names = cytidine_base_names;
        if (residue_name == "URA") base_names = uracil_base_names;
        if (residue_name == "ADE") base_names = adenine_base_names;
        if (residue_name == "GUA") base_names = guanine_base_names;

        std::string altconf_str(1, altconf);
        if (base_names.size() > 0) {
            for (gemmi::Atom& atom : residue->atoms) {
                std::string alt_name(1, atom.altloc);
                if (alt_name == altconf_str) {
                    for (const std::string& base_name : base_names) {
                        if (base_name == atom.name) {
                            base_atoms_coords.push_back(atom.pos);
                        }
                    }
                }
            }
        }
        
        if (base_atoms_coords.size() < 3) {
            throw std::runtime_error("base lsq plane has no value");
        }
    }

    float phosphate_distance_to_base_plane_gemmi(gemmi::Residue* following_res) const {
        if (!following_res) {
            throw std::runtime_error("Failed to find following phosphate");
        }

        float oop = 0.0f;
        bool found = false;

        std::string altconf_str(1, altconf);
        for (gemmi::Atom& atom : following_res->atoms) {
            std::string atm_name = atom.name;
            
            if (atm_name == " P  ") {
                std::string atom_alt(1, atom.altloc);
                if (atom_alt == altconf_str) {
                    gemmi::Position pt(atom.pos);
                    
                    if (base_atoms_coords.size() < 4) {
                        std::string m = "Failed to find base atoms. Found ";
                        m += std::to_string(base_atoms_coords.size());
                        m += " atoms. ";
                        throw std::runtime_error(m);
                    } else {
                        std::pair<double, double> oop_plus_dev =
                            lsq_plane_deviation_gemmi(base_atoms_coords, pt);
                        oop = static_cast<float>(oop_plus_dev.first);
                        found = true;
                        break;
                    }
                }
            }
        }

        if (!found) {
            throw std::runtime_error("Failed to find following phosphate");
        }
        return oop;
    }

    static std::pair<double, double> lsq_plane_deviation_gemmi(
        const std::vector<gemmi::Vec3>& points,
        const gemmi::Position& point) {

        if (points.size() < 3) {
            return std::make_pair(0.0, 0.0);
        }

        gemmi::Vec3 centroid(0, 0, 0);
        for (const auto& p : points) {
            centroid += p;
        }
        centroid /= static_cast<double>(points.size());

        gemmi::Vec3 v1 = points[1] - points[0];
        gemmi::Vec3 v2 = points[2] - points[0];
        gemmi::Vec3 normal = v1.cross(v2);
        double len = normal.length();
        if (len < 1e-10) {
            return std::make_pair(0.0, 0.0);
        }
        normal /= len;

        gemmi::Vec3 vec_to_point = static_cast<gemmi::Vec3>(point) - centroid;
        double distance = std::abs(vec_to_point.dot(normal));

        double rms = 0;
        for (const auto& p : points) {
            gemmi::Vec3 vec = p - centroid;
            double d = std::abs(vec.dot(normal));
            rms += d * d;
        }
        rms = std::sqrt(rms / static_cast<double>(points.size()));

        return std::make_pair(distance, rms);
    }
};

} // namespace coot