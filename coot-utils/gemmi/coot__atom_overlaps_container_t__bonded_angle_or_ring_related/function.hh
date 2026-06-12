#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include "geometry/protein-geometry.hh"

namespace coot {

enum bonded_atom_interaction_type {
    BAI_CLASHABLE = 0,
    BAI_BONDED    = 1,
    BAI_IGNORED   = 2
};

// ── Inline helpers ──

namespace bonded_angle_or_ring_helpers {

// Trim leading and trailing spaces from an atom name
inline std::string trim_atom_name(const std::string& name) {
    size_t first = name.find_first_not_of(" ");
    if (first == std::string::npos) return "";
    size_t last = name.find_last_not_of(" ");
    return name.substr(first, last - first + 1);
}

// are_bonded_residues: check if two residues are adjacent in sequence
inline bool are_bonded(gemmi::CRA cra_1, gemmi::CRA cra_2) {
    int seq1 = cra_1.residue->seqid.num.value;
    int seq2 = cra_2.residue->seqid.num.value;
    if (seq2 - seq1 == 1) return true;
    if (seq1 == seq2 && cra_1.residue->seqid.icode < cra_2.residue->seqid.icode) return true;
    return false;
}

// is_main_chain_p: check if atom is a main-chain atom
inline bool is_main_chain(gemmi::CRA cra) {
    std::string name = trim_atom_name(cra.atom->name);
    return name == "N" || name == "CA" || name == "C" || name == "O" || name == "CB";
}

// is_linked: check if atoms are directly connected via a Connection
inline bool is_linked(gemmi::CRA cra_1, gemmi::CRA cra_2,
                      const std::vector<gemmi::Connection>& connections) {
    for (const auto& conn : connections) {
        const auto& p1 = conn.partner1;
        const auto& p2 = conn.partner2;
        bool m1 = (p1.chain_name == cra_1.chain->name &&
                   p1.res_id.name == cra_1.residue->name &&
                   trim_atom_name(p1.atom_name) == trim_atom_name(cra_1.atom->name));
        bool m2 = (p2.chain_name == cra_2.chain->name &&
                   p2.res_id.name == cra_2.residue->name &&
                   trim_atom_name(p2.atom_name) == trim_atom_name(cra_2.atom->name));
        if ((m1 && m2) || (m2 && m1)) return true;
    }
    return false;
}

// is_ss_bonded_or_CYS_CYS_SGs
inline bool is_ss_bonded(gemmi::CRA cra_1, gemmi::CRA cra_2,
                         const std::vector<gemmi::Connection>& connections) {
    bool is_cys_1 = (cra_1.residue->name == "CYS" && trim_atom_name(cra_1.atom->name) == "SG");
    bool is_cys_2 = (cra_2.residue->name == "CYS" && trim_atom_name(cra_2.atom->name) == "SG");
    if (is_cys_1 && is_cys_2 && cra_1.residue != cra_2.residue) return true;
    return is_linked(cra_1, cra_2, connections);
}

// in_same_ring: check if two atoms in the same residue are part of the same ring
inline bool in_same_ring(gemmi::CRA cra_1, gemmi::CRA cra_2,
                         std::map<std::string, std::vector<std::vector<std::string>>>& ring_list_map) {
    std::string res_name = cra_1.residue->name;
    auto it = ring_list_map.find(res_name);
    if (it == ring_list_map.end()) return false;

    const auto& ring_atoms_vec = it->second;
    std::string atom_1 = trim_atom_name(cra_1.atom->name);
    std::string atom_2 = trim_atom_name(cra_2.atom->name);

    for (const auto& ring : ring_atoms_vec) {
        for (const auto& a : ring) {
            if (a == atom_1) {
                for (const auto& b : ring) {
                    if (b == atom_2) return true;
                }
            }
        }
    }
    return false;
}

// is_angle_related_via_link: check if two atoms in different residues could be
// angle-related via a link (disulfide / other connection)
inline bool is_angle_related_via_link(gemmi::CRA cra_1, gemmi::CRA cra_2,
                                      gemmi::Model* model_p,
                                      const std::vector<gemmi::Connection>& connections,
                                      const std::vector<std::pair<std::string, std::string>>& bonds_for_at_1,
                                      const std::vector<std::pair<std::string, std::string>>& bonds_for_at_2) {
    for (const auto& conn : connections) {
        const auto& p1 = conn.partner1;
        const auto& p2 = conn.partner2;
        std::string chain_1 = p1.chain_name;
        std::string chain_2 = p2.chain_name;
        std::string res_name_1 = p1.res_id.name;
        std::string res_name_2 = p2.res_id.name;
        std::string link_atom_1 = trim_atom_name(p1.atom_name);
        std::string link_atom_2 = trim_atom_name(p2.atom_name);

        std::string cra1_name = trim_atom_name(cra_1.atom->name);
        std::string cra2_name = trim_atom_name(cra_2.atom->name);

        // Does cra_1 match one end of the connection?
        bool cra1_matches_1 = (chain_1 == cra_1.chain->name &&
                               res_name_1 == cra_1.residue->name &&
                               link_atom_1 == cra1_name);
        bool cra1_matches_2 = (chain_2 == cra_1.chain->name &&
                               res_name_2 == cra_1.residue->name &&
                               link_atom_2 == cra1_name);

        if (cra1_matches_1) {
            for (const auto& bond : bonds_for_at_1) {
                std::string bf = trim_atom_name(bond.first);
                std::string bs = trim_atom_name(bond.second);
                if (bf == cra1_name || bs == cra1_name) {
                    std::string other = (bf == cra1_name) ? bs : bf;
                    for (const auto& bond2 : bonds_for_at_2) {
                        std::string bf2 = trim_atom_name(bond2.first);
                        std::string bs2 = trim_atom_name(bond2.second);
                        if (bf2 == cra2_name || bs2 == cra2_name) {
                            std::string other2 = (bf2 == cra2_name) ? bs2 : bf2;
                            if (other == link_atom_2 || other2 == link_atom_1) return true;
                        }
                    }
                }
            }
        }
        if (cra1_matches_2) {
            for (const auto& bond : bonds_for_at_1) {
                std::string bf = trim_atom_name(bond.first);
                std::string bs = trim_atom_name(bond.second);
                if (bf == cra1_name || bs == cra1_name) {
                    std::string other = (bf == cra1_name) ? bs : bf;
                    for (const auto& bond2 : bonds_for_at_2) {
                        std::string bf2 = trim_atom_name(bond2.first);
                        std::string bs2 = trim_atom_name(bond2.second);
                        if (bf2 == cra2_name || bs2 == cra2_name) {
                            std::string other2 = (bf2 == cra2_name) ? bs2 : bf2;
                            if (other == link_atom_1 || other2 == link_atom_2) return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

} // namespace bonded_angle_or_ring_helpers

// ── Ported free function ──

inline bonded_atom_interaction_type
bonded_angle_or_ring_related_gemmi(
    gemmi::CRA cra_1,
    gemmi::CRA cra_2,
    bool exclude_mainchain_also,
    std::map<std::string, std::vector<std::pair<std::string, std::string>>> *bonded_neighbours,
    std::map<std::string, std::vector<std::vector<std::string>>> *ring_list_map,
    gemmi::Model *model_p,
    const std::vector<gemmi::Connection> &connections,
    const protein_geometry *geom_p,
    bool ignore_water_contacts_flag)
{
    using namespace bonded_angle_or_ring_helpers;

    bonded_atom_interaction_type ait = BAI_CLASHABLE;

    gemmi::Residue *res_1 = cra_1.residue;
    gemmi::Residue *res_2 = cra_2.residue;

    if (res_1 != res_2) {
        if (are_bonded(cra_1, cra_2)) {
            if (is_main_chain(cra_1)) {
                if (is_main_chain(cra_2)) {
                    ait = BAI_BONDED;
                } else {
                    std::string res_name_2 = res_2->name;
                    if (res_name_2 == "PRO") {
                        std::string at_name_2 = trim_atom_name(cra_2.atom->name);
                        if (at_name_2 == "CD") {
                            ait = BAI_BONDED;
                        } else {
                            ait = BAI_CLASHABLE;
                        }
                    } else {
                        ait = BAI_CLASHABLE;
                    }
                }
            } else {
                if (is_main_chain(cra_2)) {
                    std::string at_name_1 = trim_atom_name(cra_1.atom->name);
                    if (at_name_1 == "CD") {
                        std::string res_name_2 = res_2->name;
                        if (res_name_2 == "PRO") {
                            ait = BAI_BONDED;
                        } else {
                            ait = BAI_CLASHABLE;
                        }
                    } else {
                        ait = BAI_CLASHABLE;
                    }
                } else {
                    ait = BAI_BONDED;
                }
            }
        } else {
            std::string res_name_1 = res_1->name;
            std::string res_name_2 = res_2->name;

            if (res_name_1 == res_name_2) {
                if (res_name_1 == "HOH") {
                    if (ignore_water_contacts_flag) {
                        ait = BAI_IGNORED;
                    } else {
                        ait = BAI_CLASHABLE;
                    }
                }
            } else {
                ait = BAI_CLASHABLE;
            }
        }
    } else {
        // ── Same residue ──
        std::string res_name = res_1->name;
        std::vector<std::pair<std::string, std::string>> bps;
        std::string atom_name_1 = cra_1.atom->name;
        std::string atom_name_2 = cra_2.atom->name;

        auto it = bonded_neighbours->find(res_name);
        if (it == bonded_neighbours->end()) {
            bps = geom_p->get_bonded_and_1_3_angles(res_name, protein_geometry::IMOL_ENC_ANY);
            (*bonded_neighbours)[res_name] = bps;
        } else {
            bps = it->second;
        }

        for (unsigned int ipr = 0; ipr < bps.size(); ipr++) {
            // Compare using trimmed names since gemmi and protein_geometry may use different padding
            std::string bp1 = trim_atom_name(bps[ipr].first);
            std::string bp2 = trim_atom_name(bps[ipr].second);
            if (trim_atom_name(atom_name_1) == bp1) {
                if (trim_atom_name(atom_name_2) == bp2) {
                    ait = BAI_BONDED;
                    break;
                }
            }
            if (trim_atom_name(atom_name_2) == bp1) {
                if (trim_atom_name(atom_name_1) == bp2) {
                    ait = BAI_BONDED;
                    break;
                }
            }
        }

        if (ait == BAI_CLASHABLE) {
            bool ringed = in_same_ring(cra_1, cra_2, *ring_list_map);
            if (ringed) ait = BAI_BONDED;
        }

        if (ait == BAI_CLASHABLE) {
            // N-terminus H check — approximate by checking if this residue
            // is the first in its chain (gemmi has no isNTerminus())
            bool is_nterminus = false;
            if (cra_1.chain && !cra_1.chain->residues.empty()) {
                if (cra_1.residue == &cra_1.chain->residues[0]) {
                    is_nterminus = true;
                }
            }

            if (is_nterminus) {
                std::string at_name_1 = trim_atom_name(cra_1.atom->name);
                std::string at_name_2 = trim_atom_name(cra_2.atom->name);
                if (at_name_1 == "H1" || at_name_1 == "H2" || at_name_1 == "H3") {
                    if (at_name_2 == "H1" || at_name_2 == "H2" || at_name_2 == "H3" ||
                        at_name_2 == "CA" || at_name_2 == "N") {
                        ait = BAI_BONDED;
                    }
                } else {
                    if (at_name_1 == "CA" || at_name_1 == "N") {
                        if (at_name_2 == "H1" || at_name_2 == "H2" || at_name_2 == "H3")
                            ait = BAI_BONDED;
                    }
                }
            }
        }
    }

    if (ait == BAI_CLASHABLE) {
        if (is_linked(cra_1, cra_2, connections) ||
            is_ss_bonded(cra_1, cra_2, connections)) {
            ait = BAI_BONDED;
        } else {
            std::vector<std::pair<std::string, std::string>> bonds_for_at_1;
            std::vector<std::pair<std::string, std::string>> bonds_for_at_2;
            std::string res_name_1 = res_1->name;
            std::string res_name_2 = res_2->name;

            auto it_1 = bonded_neighbours->find(res_name_1);
            auto it_2 = bonded_neighbours->find(res_name_2);
            if (it_1 != bonded_neighbours->end()) bonds_for_at_1 = it_1->second;
            if (it_2 != bonded_neighbours->end()) bonds_for_at_2 = it_2->second;

            if (is_angle_related_via_link(cra_1, cra_2, model_p, connections,
                                          bonds_for_at_1, bonds_for_at_2))
                ait = BAI_BONDED;
        }
    }

    return ait;
}

} // namespace coot
