#pragma once

#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include "clipper/core/cell.h"
#include "clipper/core/coords.h"
#include "coot/utils/coot-utils.hh"

namespace coot {

// Minimal ShelxIns class for the gemmi port — carries the state the
// original MMDB method reads from the enclosing object.
class ShelxIns {
public:
    bool have_cell_flag;
    clipper::Cell cell;
    std::string title;
    std::vector<std::string> sfac;
    std::vector<std::string> pre_atom_lines;
    std::vector<std::string> post_atom_lines;

    ShelxIns() : have_cell_flag(false) {}

    std::pair<int, std::string> write_ins_file_internal_gemmi(
        const gemmi::Model &model,
        const std::string &filename,
        bool mol_is_from_shelx_ins);
};

// ----------------------------------------------------------------
// Free helpers (no _gemmi port exists)
// ----------------------------------------------------------------

inline int altloc_to_part_no_impl(std::string altloc)
{
    if (altloc.empty()) return 0;
    return (altloc[0] - 'A');
}

inline int get_sfac_index_impl(const std::vector<std::string> &sfac,
                               const std::string &element_name)
{
    for (size_t i = 0; i < sfac.size(); i++) {
        if (sfac[i] == element_name) {
            return static_cast<int>(i) + 1;
        }
    }
    return -1;
}

// ----------------------------------------------------------------
// write_synthetic_pre_atom_lines_gemmi
// ----------------------------------------------------------------

inline void write_synthetic_pre_atom_lines_gemmi(
    const gemmi::Model &model,
    const clipper::Cell &cell,
    const std::vector<std::string> &sfac,
    std::ofstream &f)
{
    f << "TITL PDB->ins\n";

    // CELL line — clipper::Cell_descr uses DEGREES for angles
    double alpha_deg = cell.descr().alpha() * 180.0 / M_PI;
    double beta_deg  = cell.descr().beta()  * 180.0 / M_PI;
    double gamma_deg = cell.descr().gamma() * 180.0 / M_PI;

    f << "CELL 1.54178  ";
    f.setf(std::ios::fixed);
    f.precision(4);
    f << cell.descr().a() << " "
      << cell.descr().b() << " "
      << cell.descr().c() << " "
      << alpha_deg << " "
      << beta_deg << " "
      << gamma_deg << "\n";

    // ZERR line
    f << "ZERR 4         ";
    f.precision(4);
    f << 0.001 * cell.descr().a() << "  "
      << 0.001 * cell.descr().b() << "  "
      << 0.001 * cell.descr().c() << "  "
      << "0.0000  0.0000  0.0000\n";

    // LATT
    f << "LATT -1\n";

    // SYMM — default P 21/n symmetry
    f << "SYMM X+1/2, -Y+1/2, -Z\n";

    // SFAC
    if (!sfac.empty()) {
        f << "SFAC";
        for (size_t i = 0; i < sfac.size(); i++) {
            f << " " << sfac[i];
        }
        f << "\n";
    }

    // UNIT — count atoms per element
    std::map<std::string, int> element_counts;
    for (const gemmi::Chain &ch : model.chains) {
        for (const gemmi::Residue &res : ch.residues) {
            for (const gemmi::Atom &atom : res.atoms) {
                std::string elem = atom.element.name();
                element_counts[elem]++;
            }
        }
    }

    f << "UNIT";
    for (const auto &pair : element_counts) {
        f << " " << pair.second;
    }
    // Pad UNIT to same length as SFAC line
    int sfac_len = 4 + static_cast<int>(sfac.size() * 2) - 1;
    int unit_len = 4;
    for (const auto &pair : element_counts) {
        unit_len += 1 + static_cast<int>(std::to_string(pair.second).size());
    }
    while (unit_len < sfac_len) {
        f << " 0";
        unit_len += 2;
    }
    f << "\n";
}

// ----------------------------------------------------------------
// write_ins_file_internal_gemmi — definition
// ----------------------------------------------------------------

inline std::pair<int, std::string>
ShelxIns::write_ins_file_internal_gemmi(
    const gemmi::Model &model,
    const std::string &filename,
    bool mol_is_from_shelx_ins)
{
    int istat = 0;
    std::string message;

    double u_to_b = 8.0 * M_PI * M_PI;

    if (have_cell_flag) {
        std::ofstream f(filename.c_str());
        if (f) {

            // pre-atom lines
            if (pre_atom_lines.empty()) {
                write_synthetic_pre_atom_lines_gemmi(model, cell, sfac, f);
            } else {
                for (size_t i = 0; i < pre_atom_lines.size(); i++) {
                    f << pre_atom_lines[i] << "\n";
                }
            }

            std::string current_altloc = "";
            std::vector<std::string> afix_failure_vec;

            gemmi::UnitCell gemmi_cell(
                cell.descr().a(),
                cell.descr().b(),
                cell.descr().c(),
                cell.descr().alpha(),
                cell.descr().beta(),
                cell.descr().gamma()
            );

            // run over chains (with index tracking for resno_offset)
            int chain_idx = 0;
            for (const gemmi::Chain &chain : model.chains) {
                int resno_offset = 0;
                if (!mol_is_from_shelx_ins) {
                    resno_offset = chain_idx * 1000;
                }
                chain_idx++;

                for (const gemmi::Residue &residue : chain.residues) {
                    if (residue.atoms.empty()) continue;

                    try {
                        f << "RESI ";
                        int resno = residue.seqid.num.value;
                        int resno_out = resno + resno_offset;

                        if (resno_out < 1000) f << " ";
                        if (resno_out < 100) f << " ";
                        if (resno_out < 10) f << " ";

                        std::string chain_id = chain.name;
                        if (!chain_id.empty()) {
                            f << chain_id << ":";
                        }
                        f << resno_out << " " << residue.name << "\n";
                    } catch (const std::ios::failure &e) {
                        std::cout << "WARNING:: IOS exception caught on RESI start "
                                  << e.what() << std::endl;
                    }

                    for (const gemmi::Atom &atom : residue.atoms) {
                        try {
                            float site_occ_factor = atom.occ;
                            if (!mol_is_from_shelx_ins)
                                site_occ_factor = 11.000;

                            gemmi::Position pos(atom.pos.x, atom.pos.y, atom.pos.z);
                            gemmi::Fractional cf = gemmi_cell.fractionalize(pos);

                            std::string element_name = atom.element.name();
                            int sfac_index = get_sfac_index_impl(sfac, element_name);

                            std::string this_altloc = "";
                            if (atom.altloc != '\0' && atom.altloc != ' ') {
                                this_altloc = std::string(1, atom.altloc);
                            }

                            if (this_altloc != current_altloc) {
                                int ipart = altloc_to_part_no_impl(std::string(1, atom.altloc));
                                f << "PART    " << ipart << "\n";
                            }

                            // Isotropic B factor
                            std::string at_name = atom.name;
                            float b_factor = atom.b_iso;

                            f.setf(std::ios::fixed);
                            f.precision(7);

                            if (b_factor > 0.0) {
                                b_factor /= u_to_b;
                            }

                            f << coot::util::remove_leading_spaces(at_name)
                              << "   " << sfac_index << "  "
                              << std::setprecision(7) << std::fixed
                              << cf.x << "  " << cf.y << "   " << cf.z << " "
                              << site_occ_factor << "    "
                              << b_factor << "\n";

                            current_altloc = this_altloc;
                        } catch (const std::ios::failure &e) {
                            std::cout << "WARNING:: IOS exception caught: "
                                      << e.what() << std::endl;
                        }
                    }

                    try {
                        f << " \n"; // end of a RESI
                    } catch (const std::ios::failure &e) {
                        std::cout << "WARNING:: IOS exception caught on end of a RESI "
                                  << e.what() << std::endl;
                    }
                }
            }

            // print out first 10 AFIX misses
            int n_afix_failures = static_cast<int>(afix_failure_vec.size());
            if (n_afix_failures > 10) n_afix_failures = 10;
            for (int i = 0; i < n_afix_failures; i++)
                std::cout << afix_failure_vec[i] << std::endl;
            if (n_afix_failures > 10)
                std::cout << "WARNING:: and " << afix_failure_vec.size() - 10
                          << " more AFIX failures" << std::endl;

            if (!post_atom_lines.empty()) {
                for (size_t i = 0; i < post_atom_lines.size(); i++) {
                    try {
                        f << post_atom_lines[i] << "\n";
                    } catch (const std::ios::failure &e) {
                        std::cout << "WARNING:: IOS exception caught in post atom lines "
                                  << e.what() << std::endl;
                    }
                }
            } else {
                f << "HKLF 3\n";
                f << "END\n";
            }
        }
        f.close();
        message = "INFO:: SHELXL file ";
        message += filename;
        message += " written.";
        istat = 1;
    } else {
        std::cout << "WARNING:: no cell available... failure to write ins file."
                  << std::endl;
        message = "WARNING:: no cell available... failure to write ins file.";
    }
    return std::pair<int, std::string>(istat, message);
}

} // namespace coot