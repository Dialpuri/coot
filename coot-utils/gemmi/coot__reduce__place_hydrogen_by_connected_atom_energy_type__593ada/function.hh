#pragma once

#include <string>
#include <vector>
#include <map>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>
#include "geometry/protein-geometry.hh"

namespace coot::reduce {

/* ---- helpers (no _gemmi port exists) ---- */

inline std::string trim_str(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last  = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

inline const gemmi::Atom* find_atom_in_residue(const gemmi::Residue& res,
                                                const std::string& atom_name,
                                                const std::string& altconf) {
    std::string trimmed_name = trim_str(atom_name);
    for (const auto& atom : res.atoms) {
        if (trim_str(atom.name) != trimmed_name) continue;
        if (!altconf.empty()) {
            if (atom.altloc == altconf[0]) return &atom;
        } else {
            if (atom.altloc == ' ' || atom.altloc == '\0') return &atom;
        }
    }
    return nullptr;
}

inline gemmi::Atom*
add_hydrogen_atom_gemmi(const std::string& atom_name, const clipper::Coord_orth& pos,
                        float bf, const std::string& altconf,
                        gemmi::Residue& residue) {
    gemmi::Atom new_H;
    new_H.name = atom_name;
    new_H.element = gemmi::Element("H");
    new_H.pos = gemmi::Position(pos.x(), pos.y(), pos.z());
    new_H.occ = 1.0f;
    new_H.b_iso = bf;
    if (!altconf.empty() && altconf.length() > 0) {
        new_H.altloc = altconf[0];
    }

    std::string trimmed_atom_name = trim_str(atom_name);
    for (auto& atom : residue.atoms) {
        if (trim_str(atom.name) == trimmed_atom_name) {
            if (altconf.empty()) {
                atom.pos = new_H.pos;
                atom.occ = new_H.occ;
                atom.b_iso = new_H.b_iso;
                atom.name = atom_name;
                return &atom;
            }
            if (atom.altloc == new_H.altloc) {
                atom.pos = new_H.pos;
                atom.occ = new_H.occ;
                atom.b_iso = new_H.b_iso;
                atom.name = atom_name;
                return &atom;
            }
        }
    }

    residue.atoms.push_back(new_H);
    return &residue.atoms.back();
}

inline std::vector<std::string> get_residue_alt_confs_gemmi(const gemmi::Residue& residue) {
    std::vector<std::string> alt_confs;
    alt_confs.push_back("");
    for (const auto& atom : residue.atoms) {
        if (atom.altloc == ' ' || atom.altloc == '\0') continue;
        std::string altconf(1, atom.altloc);
        bool already = false;
        for (const auto& ac : alt_confs) {
            if (ac == altconf) { already = true; break; }
        }
        if (!already) alt_confs.push_back(altconf);
    }
    return alt_confs;
}

inline std::string get_other_H_name(unsigned int iat, const coot::dictionary_residue_restraints_t& rest) {
    std::string other_H = "";
    std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
    for (unsigned int j : neighbs) {
        if (j != iat) {
            const std::string& nb_name = rest.atom_info[j].atom_id_4c;
            if (rest.atom_info[j].is_hydrogen()) {
                other_H = nb_name;
            }
        }
    }
    return other_H;
}

inline std::map<std::string, std::vector<std::string>>
third_neighbour_map(unsigned int iat_neighb, bool exclude_hydrogen,
                    const coot::dictionary_residue_restraints_t& rest) {
    std::map<std::string, std::vector<std::string>> m;
    std::vector<unsigned int> neighbs_1 = rest.neighbours(iat_neighb, exclude_hydrogen);
    for (unsigned int i1 : neighbs_1) {
        std::string at_name_1 = rest.atom_info[i1].atom_id_4c;
        std::vector<unsigned int> neighbs_2 = rest.neighbours(i1, exclude_hydrogen);
        for (unsigned int i2 : neighbs_2) {
            if (i2 != iat_neighb) {
                std::string at_name_2 = rest.atom_info[i2].atom_id_4c;
                m[at_name_1].push_back(at_name_2);
            }
        }
    }
    return m;
}

inline std::vector<std::string> get_other_H_names(unsigned int iat, const coot::dictionary_residue_restraints_t& rest) {
    std::vector<std::string> other_Hs;
    std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
    for (unsigned int j : neighbs) {
        if (j != iat && rest.atom_info[j].is_hydrogen()) {
            other_Hs.push_back(rest.atom_info[j].atom_id_4c);
        }
    }
    return other_Hs;
}

inline double util_d2rad(double deg) {
    return deg * M_PI / 180.0;
}

/* ---- position_by_bisection ---- */

inline clipper::Coord_orth position_by_bisection(const gemmi::Atom* at_1,
                                                  const gemmi::Atom* at_2,
                                                  const gemmi::Atom* at_3,
                                                  double bl) {
    clipper::Coord_orth p1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
    clipper::Coord_orth p2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
    clipper::Coord_orth p3(at_3->pos.x, at_3->pos.y, at_3->pos.z);

    clipper::Coord_orth v12(p2.x()-p1.x(), p2.y()-p1.y(), p2.z()-p1.z());
    clipper::Coord_orth v32(p2.x()-p3.x(), p2.y()-p3.y(), p2.z()-p3.z());
    clipper::Coord_orth bisect(v12.x()+v32.x(), v12.y()+v32.y(), v12.z()+v32.z());
    double len_sq = bisect.x()*bisect.x() + bisect.y()*bisect.y() + bisect.z()*bisect.z();
    if (len_sq > 0) {
        double scale = bl / std::sqrt(len_sq);
        bisect = clipper::Coord_orth(bisect.x()*scale, bisect.y()*scale, bisect.z()*scale);
    }
    return clipper::Coord_orth(p2.x()+bisect.x(), p2.y()+bisect.y(), p2.z()+bisect.z());
}

/* ---- position_by_bond_length_angle_torsion ---- */

inline clipper::Coord_orth position_by_bond_length_angle_torsion(
    const clipper::Coord_orth& at_3,
    const clipper::Coord_orth& at_2,
    const clipper::Coord_orth& at_1,
    double bl,
    double angle_rad,
    double tor_rad) {
    clipper::Coord_orth v12(at_2.x()-at_1.x(), at_2.y()-at_1.y(), at_2.z()-at_1.z());
    clipper::Coord_orth v23(at_3.x()-at_2.x(), at_3.y()-at_2.y(), at_3.z()-at_2.z());

    double len12 = std::sqrt(v12.x()*v12.x() + v12.y()*v12.y() + v12.z()*v12.z());
    if (len12 < 1e-10) {
        return clipper::Coord_orth(at_2.x()+bl, at_2.y(), at_2.z());
    }

    double dx12 = v12.x()/len12;
    double dy12 = v12.y()/len12;
    double dz12 = v12.z()/len12;

    // Cross product: dir12 x v23
    double nx = dy12 * v23.z() - dz12 * v23.y();
    double ny = dz12 * v23.x() - dx12 * v23.z();
    double nz = dx12 * v23.y() - dy12 * v23.x();
    double norm_len = std::sqrt(nx*nx + ny*ny + nz*nz);

    if (norm_len < 1e-10) {
        // Colinear — pick arbitrary normal
        if (std::fabs(dx12) < 0.9) {
            nx = dy12 * 0 - dz12 * 0;
            ny = dz12 * 1 - dx12 * 0;
            nz = dx12 * 0 - dy12 * 1;
        } else {
            nx = dy12 * 0 - dz12 * 1;
            ny = dz12 * 0 - dx12 * 0;
            nz = dx12 * 1 - dy12 * 0;
        }
        norm_len = std::sqrt(nx*nx + ny*ny + nz*nz);
    }
    if (norm_len < 1e-10) {
        nx = 1; ny = 0; nz = 0;
    } else {
        nx /= norm_len; ny /= norm_len; nz /= norm_len;
    }

    // perpendicular = normal x dir12
    double px = ny * dz12 - nz * dy12;
    double py = nz * dx12 - nx * dz12;
    double pz = nx * dy12 - ny * dx12;

    // rel_pos = -dir12 * bl * cos(angle) + px * bl * sin(angle) * cos(tor) + nx * bl * sin(angle) * sin(tor)
    double ca = std::cos(angle_rad);
    double sa = std::sin(angle_rad);
    double ct = std::cos(tor_rad);
    double st = std::sin(tor_rad);

    double rx = -dx12 * bl * ca + px * bl * sa * ct + nx * bl * sa * st;
    double ry = -dy12 * bl * ca + py * bl * sa * ct + ny * bl * sa * st;
    double rz = -dz12 * bl * ca + pz * bl * sa * ct + nz * bl * sa * st;

    return clipper::Coord_orth(at_2.x()+rx, at_2.y()+ry, at_2.z()+rz);
}

/* ---- add_OH_H ---- */

inline void add_OH_H_gemmi(const std::string &H_name,
                            const std::string &at_name_1,
                            const std::string &at_name_2,
                            const std::string &at_name_3,
                            double bl,
                            double angle,
                            double tor_inital,
                            gemmi::Residue &residue) {
    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);
    for (size_t i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_1 = find_atom_in_residue(residue, at_name_1, alt_confs[i]);
        const gemmi::Atom* at_2 = find_atom_in_residue(residue, at_name_2, alt_confs[i]);
        const gemmi::Atom* at_3 = find_atom_in_residue(residue, at_name_3, alt_confs[i]);

        if (at_1 && at_2) {
            clipper::Coord_orth pos_1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
            clipper::Coord_orth pos_2(at_2->pos.x, at_2->pos.y, at_2->pos.z);

            clipper::Coord_orth H_pos;
            if (at_3) {
                clipper::Coord_orth pos_3(at_3->pos.x, at_3->pos.y, at_3->pos.z);
                H_pos = position_by_bond_length_angle_torsion(pos_3, pos_2, pos_1, bl,
                                                               util_d2rad(angle), util_d2rad(tor_inital));
            } else {
                H_pos = position_by_bond_length_angle_torsion(pos_1, pos_2, pos_2, bl,
                                                               util_d2rad(angle), util_d2rad(tor_inital));
            }
            float bf = static_cast<float>(at_2->b_iso);
            add_hydrogen_atom_gemmi(H_name, H_pos, bf, alt_confs[i], residue);
        }
    }
}

/* ---- add_tetrahedral_hydrogen ---- */

inline void add_tetrahedral_hydrogen_gemmi(
    const std::string& H_at_name,
    const std::string& at_central_name,
    const std::string& neighb_at_name_1,
    const std::string& neighb_at_name_2,
    const std::string& neighb_at_name_3,
    double bond_length,
    gemmi::Residue& residue) {
    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);
    for (size_t i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_central  = find_atom_in_residue(residue, at_central_name, alt_confs[i]);
        const gemmi::Atom* at_neighb_1 = find_atom_in_residue(residue, neighb_at_name_1, alt_confs[i]);
        const gemmi::Atom* at_neighb_2 = find_atom_in_residue(residue, neighb_at_name_2, alt_confs[i]);
        const gemmi::Atom* at_neighb_3 = find_atom_in_residue(residue, neighb_at_name_3, alt_confs[i]);

        if (at_central && at_neighb_1) {
            clipper::Coord_orth H_pos;
            if (at_neighb_2) {
                H_pos = position_by_bisection(at_neighb_1, at_central, at_neighb_2, bond_length);
            } else if (at_neighb_3) {
                H_pos = position_by_bisection(at_neighb_1, at_central, at_neighb_3, bond_length);
            } else {
                // Fallback: extend from central away from neighbor 1
                clipper::Coord_orth p_c(at_central->pos.x, at_central->pos.y, at_central->pos.z);
                clipper::Coord_orth p_1(at_neighb_1->pos.x, at_neighb_1->pos.y, at_neighb_1->pos.z);
                double dx = p_c.x() - p_1.x();
                double dy = p_c.y() - p_1.y();
                double dz = p_c.z() - p_1.z();
                double len = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (len > 0) {
                    double sc = bond_length / len;
                    H_pos = clipper::Coord_orth(p_c.x() + dx*sc, p_c.y() + dy*sc, p_c.z() + dz*sc);
                } else {
                    H_pos = clipper::Coord_orth(p_c.x()+bond_length, p_c.y(), p_c.z());
                }
            }
            float bf = static_cast<float>(at_central->b_iso);
            add_hydrogen_atom_gemmi(H_at_name, H_pos, bf, alt_confs[i], residue);
        }
    }
}

/* ---- add_amino_single_H ---- */

inline void add_amino_single_H_gemmi(const std::string &H_at_name,
                                      const std::string &at_name_1,
                                      const std::string &at_name_2,
                                      const std::string &at_name_3,
                                      double bl,
                                      gemmi::Residue &residue) {
    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);
    for (size_t i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_1 = find_atom_in_residue(residue, at_name_1, alt_confs[i]);
        const gemmi::Atom* at_2 = find_atom_in_residue(residue, at_name_2, alt_confs[i]);
        const gemmi::Atom* at_3 = find_atom_in_residue(residue, at_name_3, alt_confs[i]);

        if (at_1 && at_2) {
            clipper::Coord_orth H_pos;
            if (at_3) {
                H_pos = position_by_bisection(at_1, at_2, at_3, bl);
            } else {
                // Extend from at_2 away from at_1
                double dx = at_2->pos.x - at_1->pos.x;
                double dy = at_2->pos.y - at_1->pos.y;
                double dz = at_2->pos.z - at_1->pos.z;
                double len = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (len > 0) {
                    double sc = bl / len;
                    H_pos = clipper::Coord_orth(at_2->pos.x + dx*sc, at_2->pos.y + dy*sc, at_2->pos.z + dz*sc);
                } else {
                    H_pos = clipper::Coord_orth(at_2->pos.x + bl, at_2->pos.y, at_2->pos.z);
                }
            }
            float bf = static_cast<float>(at_2->b_iso);
            add_hydrogen_atom_gemmi(H_at_name, H_pos, bf, alt_confs[i], residue);
        }
    }
}

/* ---- add_2_sp3_hydrogens ---- */

inline void add_2_sp3_hydrogens_gemmi(
    const std::string &H_at_name_1,
    const std::string &H_at_name_2,
    const std::string &at_name_1,
    const std::vector<std::string> &second_neighb_vec,
    double bond_length,
    double angle_between_Hs,
    gemmi::Residue *residue_p,
    bool choose_only_farthest_position) {
    if (!residue_p) return;

    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(*residue_p);
    for (size_t i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_1 = find_atom_in_residue(*residue_p, at_name_1, alt_confs[i]);
        if (!at_1) continue;

        std::vector<const gemmi::Atom*> second_atoms;
        for (const auto& sn : second_neighb_vec) {
            const gemmi::Atom* a = find_atom_in_residue(*residue_p, sn, alt_confs[i]);
            if (a) second_atoms.push_back(a);
        }

        if (second_atoms.empty()) continue;

        float bf = static_cast<float>(at_1->b_iso);

        if (second_atoms.size() >= 2) {
            // Two H positions using bisection
            clipper::Coord_orth h1_pos = position_by_bisection(at_1, at_1, second_atoms[0], bond_length);
            clipper::Coord_orth h2_pos = position_by_bisection(at_1, at_1, second_atoms[1], bond_length);
            add_hydrogen_atom_gemmi(H_at_name_1, h1_pos, bf, alt_confs[i], *residue_p);
            add_hydrogen_atom_gemmi(H_at_name_2, h2_pos, bf, alt_confs[i], *residue_p);
        } else if (second_atoms.size() == 1) {
            clipper::Coord_orth h1_pos = position_by_bisection(at_1, at_1, second_atoms[0], bond_length);
            add_hydrogen_atom_gemmi(H_at_name_1, h1_pos, bf, alt_confs[i], *residue_p);
        }
    }
}

/* ---- main function declarations ---- */

std::vector<std::string> place_hydrogen_by_connected_atom_energy_type_gemmi(
    unsigned int iat,
    unsigned int iat_neighb,
    const coot::dictionary_residue_restraints_t& rest,
    gemmi::Residue* residue_p,
    double bl_aliph,
    double bl_arom,
    double bl_amino,
    double bl_oh,
    double bl_sh);

std::vector<std::string> place_hydrogen_by_connected_atom_energy_type_gemmi(
    const std::string& energy_type,
    unsigned int iat,
    unsigned int iat_neighb,
    const coot::dictionary_residue_restraints_t& rest,
    gemmi::Residue* residue_p,
    double bl_aliph,
    double bl_arom,
    double bl_amino,
    double bl_oh,
    double bl_sh);

} // namespace coot::reduce