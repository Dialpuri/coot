#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-map-utils.hh"
#include "clipper/clipper.h"

// Include coot-utils for float_to_string — must come AFTER clipper/clipper.h
// to avoid namespace coot shadowing std
#include "utils/coot-utils.hh"

namespace coot {

// Inline helper: build atom_spec_t from gemmi Atom+Residue+Chain
inline atom_spec_t make_atom_spec_from_gemmi(const gemmi::Chain& chain,
                                              const gemmi::Residue& res,
                                              const gemmi::Atom& atom,
                                              int model_number)
{
    atom_spec_t as;
    as.model_number = model_number;
    as.chain_id = chain.name;
    as.res_no = res.seqid.num.value;
    as.ins_code = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
    // Reconstruct original PDB atom name spacing from trimmed gemmi name
    // PDB convention: 1-char=" X  ", 2-char="XX  ", 3-char=" XXX", 4-char="XXXX"
    std::string atom_name = atom.name;
    if (atom_name.size() == 1) {
        as.atom_name = " " + atom_name + "  ";
    } else if (atom_name.size() == 2) {
        as.atom_name = atom_name + "  ";
    } else if (atom_name.size() == 3) {
        as.atom_name = " " + atom_name;
    } else {
        as.atom_name = atom_name;
    }
    as.alt_conf = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
    as.int_user_data = -1;
    as.float_user_data = -1;
    return as;
}

// Port of coot::find_water_baddies_OR to gemmi
// Finds "bad" water atoms based on B-factor, map density, and distance criteria
std::vector<atom_spec_t> find_water_baddies_OR_gemmi(
    const gemmi::Model& model,
    int model_number,
    const std::vector<gemmi::Atom*>& selected_atoms,
    float b_factor_lim,
    const clipper::Xmap<float>& xmap_in,
    float map_in_sigma,
    float outlier_sigma_level,
    float min_dist,
    float max_dist,
    short int ignore_part_occ_contact_flag,
    short int ignore_zero_occ_flag
) {
    std::vector<atom_spec_t> v;
    // Store chain+residue+atom pointers + density for each marked atom
    struct Marked {
        const gemmi::Chain* chain;
        const gemmi::Residue* residue;
        const gemmi::Atom* atom;
        float density;
    };
    std::vector<Marked> marked_for_display;

    bool this_is_marked;
    float den = 0.0f;
    short int use_b_factor_limit_test = 1;
    short int use_map_sigma_limit_test = 1;
    short int use_min_dist_test = 1;
    short int use_max_dist_test = 1;
    bool sigma_warned = false;

    if (b_factor_lim < 0.0f) use_b_factor_limit_test = 0;
    if (outlier_sigma_level < -50.0f) use_map_sigma_limit_test = 0;
    if (min_dist < 0.0f) use_min_dist_test = 0;
    if (max_dist < 0.0f) use_max_dist_test = 0;

    if (!selected_atoms.empty()) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                std::string resname = residue.name;
                if (resname == "WAT" || resname == "HOH") {
                    for (const gemmi::Atom& atom : residue.atoms) {
                        bool water_atom_is_hydrogen_atom = false;
                        if (atom.name.size() >= 2 && atom.name.compare(0, 2, " H") == 0)
                            water_atom_is_hydrogen_atom = true;
                        if (atom.name.size() >= 2 && atom.name.compare(0, 2, " D") == 0)
                            water_atom_is_hydrogen_atom = true;

                        if (water_atom_is_hydrogen_atom) continue;

                        this_is_marked = false;

                        // density check
                        if (map_in_sigma > 0.0f) {
                            clipper::Coord_orth a(atom.pos.x, atom.pos.y, atom.pos.z);
                            den = ::coot::util::density_at_point(xmap_in, a);
                            den /= map_in_sigma;
                            if (den < outlier_sigma_level && use_map_sigma_limit_test) {
                                this_is_marked = true;
                                marked_for_display.push_back({&chain, &residue, &atom, den});
                            }
                        } else {
                            if (!sigma_warned) {
                                std::cout << "Ooops! Map sigma is " << map_in_sigma << std::endl;
                                sigma_warned = true;
                            }
                        }

                        // B factor check
                        if (!this_is_marked) {
                            if (atom.b_iso > b_factor_lim && use_b_factor_limit_test) {
                                marked_for_display.push_back({&chain, &residue, &atom, den});
                            }
                        }

                        // distance check
                        if (!this_is_marked) {
                            if (ignore_part_occ_contact_flag == 0) {
                                if (ignore_zero_occ_flag == false || atom.occ > 0.01f) {
                                    double d_sqrd_min = 999999999.0;
                                    clipper::Coord_orth a(atom.pos.x, atom.pos.y, atom.pos.z);

                                    for (gemmi::Atom* other : selected_atoms) {
                                        if (&atom != other) {
                                            if (!other->element.is_hydrogen()) {
                                                clipper::Coord_orth p(other->pos.x, other->pos.y, other->pos.z);
                                                double d_sqrd = (p - a).lengthsq();
                                                if (d_sqrd < d_sqrd_min) {
                                                    d_sqrd_min = d_sqrd;
                                                }
                                            }
                                        }
                                    }

                                    bool dist_to_atoms_min_is_set = false;
                                    double dist_to_atoms_min = 99999.0;
                                    if (d_sqrd_min < 999999998.0) {
                                        dist_to_atoms_min_is_set = true;
                                        dist_to_atoms_min = std::sqrt(d_sqrd_min);
                                    }

                                    bool failed_min_dist_test = false;
                                    bool failed_max_dist_test = false;

                                    if (dist_to_atoms_min_is_set && (dist_to_atoms_min < min_dist) && use_min_dist_test)
                                        failed_min_dist_test = true;
                                    if (dist_to_atoms_min_is_set && (dist_to_atoms_min > max_dist) && use_max_dist_test)
                                        failed_max_dist_test = true;

                                    if (failed_min_dist_test || failed_max_dist_test) {
                                        marked_for_display.push_back({&chain, &residue, &atom, den});
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "marked_for_display size() " << marked_for_display.size() << std::endl;

    for (unsigned int i = 0; i < marked_for_display.size(); i++) {
        std::string s = "B fac: ";
        s += ::coot::util::float_to_string(marked_for_display[i].atom->b_iso);
        if (map_in_sigma > 0.0f) {
            s += "   ED: ";
            s += ::coot::util::float_to_string(marked_for_display[i].density);
            s += " rmsd";
        }

        atom_spec_t as = make_atom_spec_from_gemmi(*marked_for_display[i].chain,
                                                    *marked_for_display[i].residue,
                                                    *marked_for_display[i].atom,
                                                    model_number);
        as.float_user_data = marked_for_display[i].atom->occ;
        v.push_back(as);
    }

    return v;
}

} // namespace coot
