#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__get_asc_index/gemmi/function.hh"

namespace coot {

// Helper: trim leading and trailing whitespace from atom names.
// Dictionary atom names from atom_id() are 4-char padded (MMDB format),
// while gemmi stores unpadded names. This bridges the gap.
static inline std::string trim_spaces(const std::string& s) {
    auto start = s.find_first_not_of(' ');
    if (start == std::string::npos) return std::string();
    auto end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

// Gemmi port of coot::restraints_container_t::add_planes_multiatom_eigen
//
// Key changes from MMDB:
//   - mmdb::PPAtom res_selection → std::vector<gemmi::CRA> cras
//   - mmdb::PResidue SelRes → gemmi::Residue* sel_res + gemmi::Chain* sel_chain
inline int add_planes_multiatom_eigen_gemmi(
    int idr,
    const std::vector<gemmi::CRA>& cras,
    int i_no_res_atoms,
    gemmi::Residue* sel_res,
    gemmi::Chain* sel_chain,
    const coot::protein_geometry& geom)
{
    int n_plane_restr = 0;

    std::vector<std::string> altconfs = coot::util::get_residue_alt_confs_gemmi(*sel_res);
    bool all_altconfs_the_same = true;
    if (altconfs.size() > 1)
        all_altconfs_the_same = false;

    for (unsigned int ip = 0; ip < geom[idr].second.plane_restraint.size(); ip++) {
        std::map<std::string, std::vector<std::pair<int, double>>> idx_and_sigmas;

        for (int iat = 0; iat < i_no_res_atoms; iat++) {
            // gemmi stores unpadded atom names (e.g. "CA", "C1")
            std::string pdb_atom_name(cras[iat].atom->name);

            // alt_conf: "" means no alternative location, "A" means alt=A, etc.
            // gemmi stores ' ' (space) for no alt location, convert to ""
            char altloc_char = cras[iat].atom->altloc;
            std::string alt_conf;
            if (altloc_char != ' ') {
                alt_conf = std::string(1, altloc_char);
            }

            for (int irest_at = 0; irest_at < geom[idr].second.plane_restraint[ip].n_atoms(); irest_at++) {
                // Dictionary atom_id() returns 4-char padded names (MMDB format).
                // Trim to match gemmi's unpadded format.
                std::string dict_atom_name = trim_spaces(
                    geom[idr].second.plane_restraint[ip].atom_id(irest_at));

                if (pdb_atom_name == dict_atom_name) {
                    // Build atom_spec_t for get_asc_index_gemmi
                    coot::atom_spec_t spec;
                    spec.atom_name = pdb_atom_name;
                    spec.res_no    = sel_res->seqid.num.value;
                    spec.chain_id  = sel_chain->name;
                    // ins_code: gemmi uses ' ' for no insertion code → use ""
                    spec.ins_code  = "";
                    // alt_conf for get_asc_index_gemmi: it compares a->altloc against
                    // (spec.alt_conf.empty() ? '\0' : spec.alt_conf[0])
                    // gemmi stores ' ' for no alt location, so we pass " " to match
                    spec.alt_conf  = alt_conf.empty() ? " " : alt_conf;

                    int idx = coot::get_asc_index_gemmi(cras, spec);

                    if (idx >= 0) {
                        double sigma = geom[idr].second.plane_restraint[ip].dist_esd(irest_at);
                        if (sigma > 0) {
                            std::pair<int, double> idx_sigma_pair(idx, sigma);
                            if (alt_conf.empty()) {
                                if (all_altconfs_the_same) {
                                    idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
                                } else {
                                    idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
                                    for (unsigned int ialtconf = 0; ialtconf < altconfs.size(); ialtconf++) {
                                        if (!altconfs[ialtconf].empty()) {
                                            idx_and_sigmas[altconfs[ialtconf]].push_back(idx_sigma_pair);
                                        }
                                    }
                                }
                            } else {
                                idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
                            }
                        }
                    }
                }
            }
        }

        // Original calls make_fixed_flags(pos) and add_plane(..., fixed_flags)
        // For the gemmi port we count the number of plane restraints that would be added.
        for (auto it = idx_and_sigmas.begin(); it != idx_and_sigmas.end(); ++it) {
            if (it->second.size() > 3) {
                n_plane_restr++;
            }
        }
    }

    return n_plane_restr;
}

} // namespace coot