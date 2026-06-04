#pragma once

#include <gemmi/model.hpp>
#include <gemmi/symmetry.hpp>
#include <clipper/clipper.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_cell_symm/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__mol_is_anisotropic/gemmi/function.hh"

namespace coot {

// --- helper: check if two doubles are close (within 1e-6) ---
inline bool close_double_p(double a, double b) {
    return std::fabs(a - b) < 1e-6;
}

// --- helper: uppercase a string ---
inline std::string Upper(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

// --- helper: remove leading spaces from a string ---
inline std::string remove_leading_spaces(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    return s.substr(start);
}

// --- gemmi-compatible hetatom_range (stores values directly, no parent pointers) ---
struct hetatom_range_gemmi {
    std::string first_element;
    int first_seqnum;
    std::string last_element;
    int last_seqnum;
    int resno_offset;

    hetatom_range_gemmi() : first_seqnum(0), last_seqnum(0), resno_offset(0) {}
};

class ShelxIns {
public:
    std::vector<std::string> sfac;
    bool have_cell_flag = false;
    clipper::Cell cell;

    // get_atomic_contents implemented inline using gemmi primitives
    std::map<std::string, unsigned int> get_atomic_contents_gemmi(const gemmi::Structure& st) const {
        std::map<std::string, unsigned int> m;
        if (st.models.empty()) return m;
        const gemmi::Model& model = st.models[0];
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                for (const gemmi::Atom& atom : residue.atoms) {
                    std::string ele = atom.element.name();
                    if (!ele.empty()) {
                        m[ele]++;
                    }
                }
            }
        }
        return m;
    }

    void add_sfac(const std::string& element) {
        if (std::find(sfac.begin(), sfac.end(), element) == sfac.end()) {
            sfac.push_back(element);
        }
    }

    void write_synthetic_pre_atom_lines_gemmi(const gemmi::Structure& st,
                                               std::ofstream &f) {

        // TITL
        f << "TITL PDB->ins\n";

        if (have_cell_flag) {
            // Use the existing gemmi port for cell/symmetry
            std::pair<clipper::Cell, clipper::Spacegroup> cell_sg = coot::util::get_cell_symm_gemmi(st);
            const clipper::Cell& cell = cell_sg.first;
            const clipper::Spacegroup& sg = cell_sg.second;

            f << "CELL 1.54178  ";
            f << std::right << std::setprecision(4) << std::fixed
              << cell.descr().a() << " "
              << cell.descr().b() << " "
              << cell.descr().c() << " " 
              << clipper::Util::rad2d(cell.descr().alpha()) << " " 
              << clipper::Util::rad2d(cell.descr().beta())  << " " 
              << clipper::Util::rad2d(cell.descr().gamma()) << "\n";

            int n_symm = sg.num_symops();
            f << "ZERR " << n_symm << "         "
              << 0.001 * cell.descr().a() << "  " 
              << 0.001 * cell.descr().b() << "  " 
              << 0.001 * cell.descr().c() << "  ";

            // ZERR for angles
            double zerr_angle[3] = {0.05, 0.05, 0.05};
            if (close_double_p(clipper::Util::rad2d(cell.descr().alpha()),  90.0))
                zerr_angle[0] = 0.0;
            if (close_double_p(clipper::Util::rad2d(cell.descr().alpha()), 120.0))
                zerr_angle[0] = 0.0;
            if (close_double_p(clipper::Util::rad2d(cell.descr().beta()),   90.0))
                zerr_angle[1] = 0.0;
            if (close_double_p(clipper::Util::rad2d(cell.descr().beta()),  120.0))
                zerr_angle[1] = 0.0;
            if (close_double_p(clipper::Util::rad2d(cell.descr().gamma()),  90.0))
                zerr_angle[2] = 0.0;
            if (close_double_p(clipper::Util::rad2d(cell.descr().gamma()), 120.0))
                zerr_angle[2] = 0.0;
            f << zerr_angle[0] << "  "
              << zerr_angle[1] << "  "
              << zerr_angle[2] << "\n";

            // Get space group from structure
            std::string spg = st.spacegroup_hm;
            if (spg.length() > 1) {
                int latt = 1;
                char lat_char = spg[0];
                if (lat_char == 'P')
                    latt = 1;
                if (lat_char == 'I')
                    latt = 2;
                if (lat_char == 'R')
                    latt = 3;
                if (lat_char == 'F')
                    latt = 4;
                if (lat_char == 'A')
                    latt = 5;
                if (lat_char == 'B')
                    latt = 6;
                if (lat_char == 'C')
                    latt = 7;

                if (!sg.is_null()) {
                    // make latt negative if non-centrosymmetric
                    if (sg.num_inversion_symops() <= 1)
                        latt = -latt;
                }
                f << "LATT " << latt << "\n";

                // Output symmetry operations (skip the first one, isym starts from 1)
                for (int isym = 1; isym < sg.num_primitive_symops(); isym++) {
                    f << "SYMM " << Upper(sg.primitive_symop(isym).format()) << "\n";
                }
                f << "\n";
            }

            // SFAC & UNIT
            std::map<std::string, unsigned int> atomic_contents = this->get_atomic_contents_gemmi(st);
            if (!atomic_contents.empty()) {
                f << "SFAC ";
                for (const auto& kv : atomic_contents)
                    f << " " << kv.first << " ";
                f << "\n";
                f << "UNIT ";
                for (const auto& kv : atomic_contents)
                    f << kv.second * n_symm << " ";
                f << "\n";
                for (const auto& kv : atomic_contents) {
                    add_sfac(kv.first);
                }
            }
        }

        // Fixed lines
        f << "CGLS 30 -1\n";
        f << "SHEL 10 0.1\n";
        f << "FMAP 2\n";
        f << "PLAN 200 2.3\n";
        f << "LIST 6\n";
        f << "WPDB 2\n";

        bool is_aniso = coot::mol_is_anisotropic_gemmi(st);

        // --- ISOR and CONN: collect hetatom ranges ---
        std::vector<hetatom_range_gemmi> hetatom_ranges;

        if (!st.models.empty()) {
            const gemmi::Model& model = st.models[0];
            for (unsigned int ichain = 0; ichain < model.chains.size(); ichain++) {
                int resno_offset = static_cast<int>(ichain) * 1000;
                const gemmi::Chain& chain = model.chains[ichain];

                std::string running_hetatm_elem;
                int running_hetatm_seqnum = 0;
                hetatom_range_gemmi current_range;

                for (const gemmi::Residue& residue : chain.residues) {
                    for (const gemmi::Atom& atom : residue.atoms) {
                        // Skip hydrogen atoms (replacement for MMDB isTer())
                        if (!atom.element.is_hydrogen()) {
                            if (residue.het_flag != ' ') {
                                running_hetatm_elem = atom.element.name();
                                running_hetatm_seqnum = residue.seqid.num.value;
                                if (current_range.first_element.empty()) {
                                    current_range.first_element = running_hetatm_elem;
                                    current_range.first_seqnum = running_hetatm_seqnum;
                                    current_range.resno_offset = resno_offset;
                                }
                            } else {
                                if (!current_range.first_element.empty()) {
                                    current_range.last_element = running_hetatm_elem;
                                    current_range.last_seqnum = running_hetatm_seqnum;
                                    hetatom_ranges.push_back(current_range);
                                    current_range = hetatom_range_gemmi();
                                }
                            }
                        }
                    }
                }

                // End of chain: if we have an open range, close it
                if (!running_hetatm_elem.empty() && !current_range.first_element.empty()) {
                    if (current_range.last_element.empty()) {
                        current_range.last_element = running_hetatm_elem;
                        current_range.last_seqnum = running_hetatm_seqnum;
                        hetatom_ranges.push_back(current_range);
                    }
                }
            }
        }

        // Write ISOR for anisotropic case
        if (is_aniso) {
            for (const auto& hr : hetatom_ranges) {
                f << "ISOR 0.1 "
                  << remove_leading_spaces(hr.first_element) << "_"
                  << hr.first_seqnum + hr.resno_offset << " > " 
                  << remove_leading_spaces(hr.last_element) << "_" 
                  << hr.last_seqnum + hr.resno_offset << "\n";
            }
        }

        // Write CONN
        for (const auto& hr : hetatom_ranges) {
            f << "CONN 0 "
              << remove_leading_spaces(hr.first_element) << "_"
              << hr.first_seqnum + hr.resno_offset << " > " 
              << remove_leading_spaces(hr.last_element) << "_" 
              << hr.last_seqnum + hr.resno_offset << "\n";
        }

        // DELU and RIGU for anisotropic
        if (is_aniso) {
            std::map<std::string, unsigned int> atomic_contents = this->get_atomic_contents_gemmi(st);
            if (!atomic_contents.empty()) {
                f << "DELU ";
                for (const auto& kv : atomic_contents)
                    f << " %" << kv.first << "_* ";
                f << "\n";
            }
            f << "RIGU\n";
        }

        f << "BUMP\n";
    }
};

} // namespace coot