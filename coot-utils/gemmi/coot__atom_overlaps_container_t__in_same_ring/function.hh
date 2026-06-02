#pragma once
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <gemmi/model.hpp>
#include "coot-utils/atom-overlaps.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

// gemmi port of atom_overlaps_container_t::in_same_ring
// Free function because atom_overlaps_container_t members are private.
// Takes gemmi::CRA instead of mmdb::Atom* for parent context access.
inline bool in_same_ring_gemmi(
    const protein_geometry* geom,
    gemmi::CRA cra_1,
    gemmi::CRA cra_2,
    std::map<std::string, std::vector<std::vector<std::string>>> &ring_list_map
) {
    bool same = false;

    if (cra_1.residue == nullptr || cra_2.residue == nullptr) {
        return false;
    }

    if (cra_1.residue == cra_2.residue) {
        // Trim atom names (gemmi stores them padded to 4 chars)
        std::string at_name_1 = cra_1.atom->name;
        at_name_1.erase(at_name_1.find_last_not_of(" \t\r\n") + 1);

        std::string at_name_2 = cra_2.atom->name;
        at_name_2.erase(at_name_2.find_last_not_of(" \t\r\n") + 1);

        try {
            // residue->name is already std::string in gemmi
            std::string res_name = cra_1.residue->name;

            // Get dictionary from geometry
            auto p = geom->get_monomer_restraints(res_name, protein_geometry::IMOL_ENC_ANY);
            if (!p.first) {
                return false;
            }
            const dictionary_residue_restraints_t &dict = p.second;

            std::map<std::string, std::vector<std::vector<std::string>>>::const_iterator it;
            it = ring_list_map.find(res_name);
            if (it != ring_list_map.end()) {
                same = dict.in_same_ring(at_name_1, at_name_2, it->second);
            } else {
                std::vector<std::vector<std::string>> ring_list;
                if (res_name == "HIS") {
                    ring_list = {{"ND1", "CD2", "CE1", "NE2"}};
                } else if (res_name == "PHE" || res_name == "TYR") {
                    ring_list = {{"CD1", "CE1", "CZ", "CE2", "CD2", "CG"}};
                } else if (res_name == "TRP") {
                    ring_list = {
                        {"CD1", "NE1", "CE2", "CZ2", "CH2", "CE3", "CZ3", "CD2", "CG"},
                        {"CD1", "NE1", "CE2", "CZ2", "CH2", "CE3", "CZ3", "CD2"}
                    };
                } else if (res_name == "PRO") {
                    ring_list = {{"N", "CA", "CB", "CG", "CD"}};
                } else {
                    ring_list = dict.get_ligand_ring_list();
                }
                ring_list_map[res_name] = ring_list;
                same = dict.in_same_ring(at_name_1, at_name_2, ring_list);
            }
        }
        catch (const std::out_of_range &ex) {
            std::cout << "OOpps " << ex.what() << std::endl;
        }
    }
    return same;
}

} // namespace coot