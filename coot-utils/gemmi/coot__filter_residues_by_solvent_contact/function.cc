#include "function.hh"
#include <cmath>

namespace coot {

std::vector<const gemmi::Residue*>
filter_residues_by_solvent_contact_gemmi(const gemmi::Residue& res_ref,
                                         const gemmi::Structure& st,
                                         const std::vector<const gemmi::Residue*>& residues,
                                         double water_dist_max) {
    std::vector<const gemmi::Residue*> v;
    const double max_dist_sq = water_dist_max * water_dist_max;

    // Get non-H/C atoms from reference residue
    std::vector<gemmi::Position> ligand_non_hc_atoms;
    for (const gemmi::Atom& atom : res_ref.atoms) {
        std::string ele = atom.element.name();
        if (ele != "C" && ele != "H") {
            ligand_non_hc_atoms.push_back(atom.pos);
        }
    }

    for (const gemmi::Residue* res : residues) {
        std::string res_name = res->name;
        if (res_name != "HOH") {
            v.push_back(res);
        } else {
            bool added = false;
            for (const gemmi::Atom& atom : res->atoms) {
                if (added) break;
                gemmi::Position pt = atom.pos;
                for (const gemmi::Position& lig_pt : ligand_non_hc_atoms) {
                    double dist_sq = (lig_pt - pt).length_sq();
                    if (dist_sq < max_dist_sq) {
                        added = true;
                        v.push_back(res);
                        break;
                    }
                }
            }
        }
    }

    return v;
}

} // namespace coot