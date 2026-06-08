#include "function.hh"

namespace coot {

bool restraints_container_t::add_torsion_internal_gemmi(
    const dict_torsion_restraint_t &torsion_restraint,
    const std::vector<gemmi::Atom*> &res_atoms,
    const std::map<gemmi::Atom*, int> &atom_index_map,
    const double &torsion_restraint_weight)
{
    bool status = false;

    // Helper: pad atom name to 4 chars with MMDB/PDB convention
    // (left-pad with one space, then right-pad to 4 total)
    // "N"  → " N  ", "CA" → " CA ", "CB" → " CB "
    auto pad4 = [](const std::string& s) {
        std::string p = " " + s;
        while (p.size() < 4) p += ' ';
        return p;
    };

    // Joel Bard fix: Don't add torsion restraints for torsion that
    // have either s.d. or period 0
    if (torsion_restraint.periodicity() > 0) {
        if (torsion_restraint.esd() > 0.000001) {

            int n = static_cast<int>(res_atoms.size());

            for (int iat = 0; iat < n; iat++) {
                std::string pdb_atom_name1(pad4(res_atoms[iat]->name));

                if (pdb_atom_name1 == torsion_restraint.atom_id_1_4c()) {
                    for (int iat2 = 0; iat2 < n; iat2++) {
                        std::string pdb_atom_name2(pad4(res_atoms[iat2]->name));
                        if (pdb_atom_name2 == torsion_restraint.atom_id_2_4c()) {
                            for (int iat3 = 0; iat3 < n; iat3++) {
                                std::string pdb_atom_name3(pad4(res_atoms[iat3]->name));
                                if (pdb_atom_name3 == torsion_restraint.atom_id_3_4c()) {
                                    for (int iat4 = 0; iat4 < n; iat4++) {
                                        std::string pdb_atom_name4(pad4(res_atoms[iat4]->name));
                                        if (pdb_atom_name4 == torsion_restraint.atom_id_4_4c()) {

                                            // kill off weird dictionary torsions here
                                            if (pdb_atom_name1 == " O  ")
                                                if (pdb_atom_name2 == " C  ")
                                                    if (pdb_atom_name3 == " CA ")
                                                        continue;
                                            if (pdb_atom_name1 == " CB ")
                                                if (pdb_atom_name2 == " CA  ")
                                                    if (pdb_atom_name3 == " N ")
                                                        if (pdb_atom_name4 == " H ")
                                                            continue;

                                            // Look up atom indices from the map
                                            int index1 = 0, index2 = 0, index3 = 0, index4 = 0;
                                            {
                                                auto it = atom_index_map.find(res_atoms[iat]);
                                                if (it != atom_index_map.end()) index1 = it->second;
                                            }
                                            {
                                                auto it = atom_index_map.find(res_atoms[iat2]);
                                                if (it != atom_index_map.end()) index2 = it->second;
                                            }
                                            {
                                                auto it = atom_index_map.find(res_atoms[iat3]);
                                                if (it != atom_index_map.end()) index3 = it->second;
                                            }
                                            {
                                                auto it = atom_index_map.find(res_atoms[iat4]);
                                                if (it != atom_index_map.end()) index4 = it->second;
                                            }

                                            double torsion_angle = torsion_restraint.angle();
                                            if (torsion_angle < 0)
                                                torsion_angle += 360;
                                            if (torsion_angle > 360)
                                                torsion_angle -= 360;

                                            // gemmi altloc: ' ' (space) = no alt loc, matches MMDB ""
                                            char c1 = res_atoms[iat]->altloc;
                                            char c2 = res_atoms[iat2]->altloc;
                                            char c3 = res_atoms[iat3]->altloc;
                                            char c4 = res_atoms[iat4]->altloc;

                                            std::string a1 = (c1 == ' ' ? "" : std::string(1, c1));
                                            std::string a2 = (c2 == ' ' ? "" : std::string(1, c2));
                                            std::string a3 = (c3 == ' ' ? "" : std::string(1, c3));
                                            std::string a4 = (c4 == ' ' ? "" : std::string(1, c4));

                                            bool alt_confs_match = false;
                                            if (a1 == "" || a1 == a2)
                                                if (a2 == "" || a2 == a3)
                                                    if (a3 == "" || a3 == a4 || a4 == "")
                                                        alt_confs_match = true;

                                            if (alt_confs_match) {
                                                std::vector<bool> fixed_flags = make_fixed_flags(index1, index2, index3, index4);
                                                add(0, index1, index2, index3, index4,
                                                    fixed_flags,
                                                    torsion_angle,
                                                    torsion_restraint.esd(),
                                                    torsion_restraint_weight,
                                                    torsion_restraint.periodicity());
                                                status = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return status;
}

}  // namespace coot
