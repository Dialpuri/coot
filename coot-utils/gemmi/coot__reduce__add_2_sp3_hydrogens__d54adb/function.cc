#include "function.hh"

namespace coot { namespace reduce {

void add_2_sp3_hydrogens_gemmi(const std::string &H_at_name_1,
                               const std::string &H_at_name_2,
                               const std::string &at_name_1,
                               const std::string &at_name_2,
                               const std::string &at_name_3,
                               double bond_length,
                               double angle_between_Hs,
                               gemmi::Residue &residue,
                               bool choose_only_farthest_position) {
    
    std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
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
            std::pair<clipper::Coord_orth, clipper::Coord_orth> Hs =
                position_pair_by_bisection_gemmi(at_1, at_2, at_3,
                                                 bond_length,
                                                 clipper::Util::d2rad(angle_between_Hs));
            
            float bf = at_2->b_iso;
            
            if (!choose_only_farthest_position) {
                add_hydrogen_atom_gemmi(H_at_name_1, Hs.first, bf, alt_conf, residue);
                add_hydrogen_atom_gemmi(H_at_name_2, Hs.second, bf, alt_conf, residue);
            } else {
                clipper::Coord_orth at_pos_1 = co_gemmi(at_1);
                clipper::Coord_orth at_pos_3 = co_gemmi(at_3);
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

}} // namespace coot::reduce