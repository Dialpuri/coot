#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <map>
#include "geometry/protein-geometry.hh"

namespace coot {

// Lightweight restraints_container_t for gemmi port.
// We do NOT include simple-restraint.hh (which defines its own
// coot::restraints_container_t with a different layout) to avoid
// type conflicts. This lightweight struct provides what add_torsions_gemmi needs.
struct restraints_container_t {
    std::vector<std::vector<int>> restraints_vec;

    int add_torsions_gemmi(
        int idr,
        const std::vector<gemmi::Atom*> &res_atoms,
        const std::map<gemmi::Atom*, int> &atom_index_map,
        const coot::protein_geometry &geom,
        const double &torsion_restraints_weight)
    {
        int n_torsion_restr = 0;
        const std::vector<dict_torsion_restraint_t> &torsion_restraints = geom[idr].second.torsion_restraint;

        // Helper: pad atom name to 4 chars with PDB convention (left-pad space, right-pad to 4)
        auto pad4 = [](const std::string& s) {
            std::string p = " " + s;
            while (p.size() < 4) p += ' ';
            return p;
        };

        for (unsigned int ib = 0; ib < torsion_restraints.size(); ib++) {
            const dict_torsion_restraint_t &tr = torsion_restraints[ib];

            // Joel Bard fix: skip torsions with period 0 or esd ~0
            if (tr.periodicity() <= 0) continue;
            if (tr.esd() <= 0.000001) continue;

            int n = static_cast<int>(res_atoms.size());

            for (int iat = 0; iat < n; iat++) {
                std::string pdb_atom_name1 = pad4(res_atoms[iat]->name);
                if (pdb_atom_name1 != tr.atom_id_1_4c()) continue;

                for (int iat2 = 0; iat2 < n; iat2++) {
                    std::string pdb_atom_name2 = pad4(res_atoms[iat2]->name);
                    if (pdb_atom_name2 != tr.atom_id_2_4c()) continue;

                    for (int iat3 = 0; iat3 < n; iat3++) {
                        std::string pdb_atom_name3 = pad4(res_atoms[iat3]->name);
                        if (pdb_atom_name3 != tr.atom_id_3_4c()) continue;

                        for (int iat4 = 0; iat4 < n; iat4++) {
                            std::string pdb_atom_name4 = pad4(res_atoms[iat4]->name);
                            if (pdb_atom_name4 != tr.atom_id_4_4c()) continue;

                            // Kill off weird dictionary torsions
                            if (pdb_atom_name1 == " O  " && pdb_atom_name2 == " C  " && pdb_atom_name3 == " CA ")
                                continue;
                            if (pdb_atom_name1 == " CB " && pdb_atom_name2 == " CA  " &&
                                pdb_atom_name3 == " N  " && pdb_atom_name4 == " H  ")
                                continue;

                            // Look up atom indices from the map
                            int index1 = 0, index2 = 0, index3 = 0, index4 = 0;
                            { auto it = atom_index_map.find(res_atoms[iat]);   if (it != atom_index_map.end()) index1 = it->second; }
                            { auto it = atom_index_map.find(res_atoms[iat2]);  if (it != atom_index_map.end()) index2 = it->second; }
                            { auto it = atom_index_map.find(res_atoms[iat3]);  if (it != atom_index_map.end()) index3 = it->second; }
                            { auto it = atom_index_map.find(res_atoms[iat4]);  if (it != atom_index_map.end()) index4 = it->second; }

                            double torsion_angle = tr.angle();
                            if (torsion_angle < 0)  torsion_angle += 360;
                            if (torsion_angle > 360) torsion_angle -= 360;

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
                                // Type 0 = torsion in the generated port convention
                                restraints_vec.push_back({0, index1, index2, index3, index4});
                                n_torsion_restr++;
                            }
                        }
                    }
                }
            }
        }

        return n_torsion_restr;
    }
};

}  // namespace coot
