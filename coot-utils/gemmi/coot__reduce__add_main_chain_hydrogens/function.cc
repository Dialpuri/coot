#include "function.hh"
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_pair_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_tetrahedron/gemmi/function.hh"

namespace coot::reduce {

// Helper: find atom by unpadded name and alt_conf
// gemmi stores atom names WITHOUT padding (e.g. "CA" not " CA ")
static const gemmi::Atom* find_atom_in_residue(const gemmi::Residue& residue,
                                                const std::string& atom_name,
                                                const std::string& altconf) {
    for (const auto& atom : residue.atoms) {
        if (atom.name == atom_name) {
            if (altconf.empty() || altconf == " ") {
                if (atom.altloc == '\0' || atom.altloc == ' ')
                    return &atom;
            } else {
                if (atom.altloc == static_cast<char>(altconf[0]))
                    return &atom;
            }
        }
    }
    return nullptr;
}

// Implement add_main_chain_HA inline — matches the ported implementation
// Uses C, N, CB as tetrahedron neighbors (not N, C, N)
static void add_main_chain_HA_impl(gemmi::Residue& residue) {
    double bl = 0.97;

    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);

    if (alt_confs.empty()) {
        alt_confs.push_back("");
    }

    for (unsigned int i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_ca  = find_atom_in_residue(residue, "CA", alt_confs[i]);
        const gemmi::Atom* at_n1  = find_atom_in_residue(residue, "C",  alt_confs[i]);
        const gemmi::Atom* at_n2  = find_atom_in_residue(residue, "N",  alt_confs[i]);
        const gemmi::Atom* at_n3  = find_atom_in_residue(residue, "CB", alt_confs[i]);

        if (at_ca && at_n1 && at_n2 && at_n3) {
            clipper::Coord_orth pos = position_by_tetrahedron_gemmi(at_ca, at_n1, at_n2, at_n3, bl);
            float bf = static_cast<float>(at_ca->b_iso);
            add_hydrogen_atom_gemmi(" HA ", pos, bf, alt_confs[i], residue);
        }
    }
}

// Implement add_2_sp3_hydrogens (3-atom-name variant) inline
static void add_2_sp3_hydrogens_impl(
    const std::string& H_at_name_1,
    const std::string& H_at_name_2,
    const std::string& at_name_1,
    const std::string& at_name_2,
    const std::string& at_name_3,
    double bond_length,
    double angle_between_Hs,
    gemmi::Residue& residue,
    bool choose_only_farthest_position) {

    std::vector<std::string> alt_confs = get_residue_alt_confs_gemmi(residue);
    if (alt_confs.empty()) alt_confs.push_back(" ");

    for (const auto& alt_conf : alt_confs) {
        const gemmi::Atom* at_1 = nullptr;
        const gemmi::Atom* at_2 = nullptr;
        const gemmi::Atom* at_3 = nullptr;

        for (const auto& atom : residue.atoms) {
            if (atom.name == at_name_1) at_1 = &atom;
            if (atom.name == at_name_2) at_2 = &atom;
            if (atom.name == at_name_3) at_3 = &atom;
        }

        if (at_1 && at_2 && at_3) {
            // Convert angle to radians before passing to position_pair_by_bisection
            std::pair<clipper::Coord_orth, clipper::Coord_orth> Hs =
                position_pair_by_bisection_gemmi(at_1, at_2, at_3,
                                                 bond_length,
                                                 clipper::Util::d2rad(angle_between_Hs));

            float bf = at_2->b_iso;

            if (!choose_only_farthest_position) {
                add_hydrogen_atom_gemmi(H_at_name_1, Hs.first,  bf, alt_conf, residue);
                add_hydrogen_atom_gemmi(H_at_name_2, Hs.second, bf, alt_conf, residue);
            } else {
                clipper::Coord_orth at_pos_1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
                clipper::Coord_orth at_pos_3(at_3->pos.x, at_3->pos.y, at_3->pos.z);
                clipper::Coord_orth mp(0.5 * (at_pos_1 + at_pos_3));
                double d1 = Hs.first.lengthsq();
                double d2 = Hs.second.lengthsq();
                if (d1 > d2)
                    add_hydrogen_atom_gemmi(H_at_name_1, Hs.first, bf, alt_conf, residue);
                else
                    add_hydrogen_atom_gemmi(H_at_name_1, Hs.second, bf, alt_conf, residue);
            }
        }
    }
}

void add_main_chain_hydrogens_gemmi(gemmi::Residue& residue, gemmi::Residue* residue_prev, bool is_gly) {
    if (is_gly) {
        double bl = 0.97;
        add_2_sp3_hydrogens_impl("HA2", "HA3", "N", "CA", "C", bl, 107, residue, false);
        add_main_chain_H_gemmi(residue, residue_prev);
    } else {
        add_main_chain_HA_impl(residue);
        std::string residue_name = residue.name;
        if (is_standard_amino_acid_name(residue_name))
            if (residue_name != "PRO")
                add_main_chain_H_gemmi(residue, residue_prev);
    }
}

} // namespace coot::reduce