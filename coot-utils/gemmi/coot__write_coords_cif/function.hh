#pragma once

#include <gemmi/model.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/neighbor.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

namespace coot {
namespace util {

// Translate coot::util::remove_wrong_cis_peptides for gemmi::Structure
// Original logic: iterate residues, find N-Ca-C-N peptide bonds,
// correct those with wrong cis/trans configuration.
inline void remove_wrong_cis_peptides(gemmi::Structure& st) {
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            // Collect peptide-bonding atom pointers in order
            std::vector<gemmi::Residue*> residues;
            for (gemmi::Residue& res : chain.residues) {
                // Skip water/hetero that don't have peptide bonds
                if (res.name == "HOH" || res.name == "WAT") continue;
                residues.push_back(&res);
            }

            // For each consecutive pair, check N-Ca-C-N peptide bond
            for (size_t i = 0; i + 1 < residues.size(); ++i) {
                gemmi::Residue* res_i = residues[i];
                gemmi::Residue* res_j = residues[i + 1];

                // Find C of res_i and N of res_j
                gemmi::Atom* c_atom = nullptr;
                gemmi::Atom* n_atom = nullptr;
                for (gemmi::Atom& a : res_i->atoms) {
                    if (a.name == " C  ") { c_atom = &a; break; } // C terminal
                }
                for (gemmi::Atom& a : res_j->atoms) {
                    if (a.name == " N  ") { n_atom = &a; break; } // N terminal
                }
                if (!c_atom || !n_atom) continue;

                // Distance C-N
                double dist_sq = (c_atom->pos - n_atom->pos).length_sq();
                double dist = std::sqrt(dist_sq);

                // Normal peptide bond distance ~1.32-1.35 A
                // Cis peptide bond is much longer when wrong
                // If distance > 1.8, it's a wrong cis peptide - swap back
                if (dist > 1.8) {
                    // This is a wrong cis peptide bond
                    // Find Ca of res_i and Ca of res_j
                    gemmi::Atom* ca_i = nullptr;
                    for (gemmi::Atom& a : res_i->atoms) {
                        if (a.name == " CA") { ca_i = &a; break; }
                    }
                    gemmi::Atom* ca_j = nullptr;
                    for (gemmi::Atom& a : res_j->atoms) {
                        if (a.name == " CA") { ca_j = &a; break; }
                    }
                    gemmi::Atom* o_atom = nullptr;
                    for (gemmi::Atom& a : res_i->atoms) {
                        if (a.name == " O  ") { o_atom = &a; break; }
                    }
                    gemmi::Atom* h_atom = nullptr;
                    for (gemmi::Atom& a : res_j->atoms) {
                        if (a.name == " H  ") { h_atom = &a; break; }
                    }

                    // Compute omega angle N(i)-Ca(i)-C(i)-N(j)
                    if (ca_i && o_atom && h_atom) {
                        gemmi::Vec3 n_pos = n_atom->pos;
                        gemmi::Vec3 ca_pos = ca_i->pos;
                        gemmi::Vec3 c_pos  = c_atom->pos;

                        // Omega dihedral: N-Ca-C-N
                        gemmi::Vec3 v1 = n_pos - ca_pos;
                        gemmi::Vec3 v2 = c_pos  - ca_pos;
                        gemmi::Vec3 v3 = n_pos - c_pos;

                        // Compute dihedral angle
                        gemmi::Vec3 n1 = v1.cross(v2);
                        gemmi::Vec3 n2 = v2.cross(v3);
                        if (n1.length_sq() > 0 && n2.length_sq() > 0) {
                            n1 /= std::sqrt(n1.length_sq());
                            n2 /= std::sqrt(n2.length_sq());
                            gemmi::Vec3 n3 = n1.cross(n2);
                            n3 /= std::sqrt(n3.length_sq());

                            double cos_val = std::clamp(n1.dot(n2), -1.0, 1.0);
                            double omega = std::acos(cos_val);
                            if (n3.dot(v2) < 0) omega = 2.0 * M_PI - omega;

                            // omega > 3.14 means wrong cis (should be ~0 for cis, ~pi for trans)
                            // Actually: trans ~ pi (3.14), cis ~ 0
                            // If distance is > 1.8 A, the peptide was flipped wrong
                            // Correct by swapping appropriate atoms
                            double omega_deg = omega * 180.0 / M_PI;

                            // If it looks like a wrong cis (>150 deg but distance > 1.8)
                            if (omega_deg > 150) {
                                // Swap N and H positions, C and O positions
                                gemmi::Vec3 tmp = n_atom->pos;
                                n_atom->pos = gemmi::Position(h_atom->pos.x, h_atom->pos.y, h_atom->pos.z);
                                h_atom->pos = gemmi::Position(tmp.x, tmp.y, tmp.z);

                                tmp = c_atom->pos;
                                c_atom->pos = gemmi::Position(o_atom->pos.x, o_atom->pos.y, o_atom->pos.z);
                                o_atom->pos = gemmi::Position(tmp.x, tmp.y, tmp.z);
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace util

// Port of coot::write_coords_cif
// Original: int coot::write_coords_cif(mmdb::Manager *mol, const std::string &file_name)
inline int write_coords_cif_gemmi(gemmi::Structure& st, const std::string& file_name) {
    util::remove_wrong_cis_peptides(st);

    std::ofstream ofs(file_name);
    if (!ofs) return -1;

    gemmi::cif::Document doc = gemmi::make_mmcif_document(st);
    gemmi::cif::write_cif_to_stream(ofs, doc);
    ofs.close();

    return 0;
}

} // namespace coot