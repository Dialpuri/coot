#include "function.hh"

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <gemmi/mmread.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/cell.h>
#include <clipper/core/spacegroup.h>
#include <clipper/contrib/edcalc.h>
#include <coot-utils/coot-density-stats.hh>
#include <coot/geometry/residue-and-atom-specs.hh>
#include <coot-utils/coot-map-utils.hh>
#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cctype>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

// Trim whitespace from both ends
static std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

// Compare unpadded atom name (handles both PDB-padded and CIF-unpadded)
static bool atom_name_is(const gemmi::Atom &at, const char *name) {
    return trim(at.name) == name;
}

// Check if residue is a het residue (het_flag == 'H' means HETATM)
static bool is_het_residue_gemmi(const gemmi::Residue *residue_p) {
    return residue_p->het_flag == 'H';
}

// Calculate electron density map from atoms in a gemmi::Model
static clipper::Xmap<float> calc_atom_map_from_model(
    const gemmi::Model &model,
    const clipper::Cell &cell,
    const clipper::Spacegroup &space_group,
    const clipper::Grid_sampling &sampling) {

    clipper::Xmap<float> xmap;
    xmap.init(space_group, cell, sampling);

    std::vector<clipper::Atom> clipper_atoms;
    float rescale_b_u = 1.0f / (8.0f * M_PI * M_PI);

    for (const gemmi::Chain &chain : model.chains) {
        for (const gemmi::Residue &residue : chain.residues) {
            for (const gemmi::Atom &at : residue.atoms) {
                clipper::Coord_orth pt(at.pos.x, at.pos.y, at.pos.z);
                std::string ele = at.element.name();
                clipper::Atom cat;
                cat.set_element(ele);
                cat.set_coord_orth(pt);
                float u_iso = at.b_iso * rescale_b_u;
                if (u_iso < 0.1f) u_iso = 0.1f;
                cat.set_u_iso(u_iso);
                cat.set_occupancy(at.occ);
                clipper_atoms.push_back(cat);
            }
        }
    }

    try {
        clipper::Atom_list al(clipper_atoms);
        clipper::EDcalc_iso<float> ed_calc;
        ed_calc(xmap, al);

        // Clean NaN/inf
        clipper::Xmap<float>::Map_reference_index inx;
        for (inx = xmap.first(); !inx.last(); inx.next()) {
            if (std::isnan(xmap[inx]) || std::isinf(xmap[inx])) {
                xmap[inx] = 0.0f;
            }
        }
    } catch (const clipper::Message_generic &e) {
        std::cout << "ERROR:: calc_atom_map_from_model clipper problem" << std::endl;
        std::cout << e.text() << std::endl;
    }
    return xmap;
}

// Local residue_run_t for gemmi (uses const gemmi::Residue*)
class residue_run_t_gemmi {
public:
    unsigned int idx_mid;
    unsigned int n_residues_per_blob;
    std::vector<const gemmi::Residue *> residues;

    residue_run_t_gemmi() : idx_mid(0), n_residues_per_blob(0) {}

    explicit residue_run_t_gemmi(unsigned int n_residues_per_blob_in)
        : n_residues_per_blob(n_residues_per_blob_in) {
        idx_mid = n_residues_per_blob / 2;
    }

    explicit residue_run_t_gemmi(const std::vector<const gemmi::Residue *> &rr_in)
        : residues(rr_in) {
        n_residues_per_blob = rr_in.size();
        idx_mid = residues.size() / 2;
    }

    void add(const std::vector<const gemmi::Residue *> &rv) {
        residues.insert(residues.begin(), rv.begin(), rv.end());
    }

    const gemmi::Residue *residue_mid() const {
        if (residues.size() >= n_residues_per_blob) {
            return residues[idx_mid];
        } else {
            std::cout << "ERROR:: indexing residues_mid " << residues.size() << std::endl;
            return nullptr;
        }
    }
};

std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
map_to_model_correlation_stats_per_residue_run_gemmi(
    const gemmi::Model &model,
    const std::string &chain_id,
    const clipper::Xmap<float> &xmap,
    unsigned int n_residues_per_blob,
    bool exclude_CON,
    float atom_mask_radius,
    float NOC_mask_radius) {

    std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> res_map_all_atom;
    std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> res_map_side_chain;

    // Collect ALL residues from the matching chain (including filtered ones)
    // This matches MMDB behavior: iterate over chain positions, not just valid residues
    std::vector<const gemmi::Residue *> all_chain_residues;
    for (const gemmi::Chain &chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (const gemmi::Residue &res : chain.residues) {
            all_chain_residues.push_back(&res);
        }
    }

    if (all_chain_residues.empty()) {
        return std::make_pair(res_map_all_atom, res_map_side_chain);
    }

    // Build residue runs using sliding window over chain positions (matching MMDB original)
    // MMDB: for (int i_res=0; i_res<idx_last_residue; i_res++) where idx_last_residue = n_residues - n_residues_per_blob
    std::vector<residue_run_t_gemmi> residue_runs;
    int n_res = static_cast<int>(all_chain_residues.size());
    int idx_last_residue = n_res - static_cast<int>(n_residues_per_blob);
    std::vector<const gemmi::Residue *> residue_vec;

    for (int i_res = 0; i_res < idx_last_residue; i_res++) {
        for (unsigned int i_run = 0; i_run < n_residues_per_blob; i_run++) {
            const gemmi::Residue *residue_p = all_chain_residues[i_res + static_cast<int>(i_run)];
            std::string res_name = residue_p->name;
            if (res_name == "HOH" ||
                res_name == "A" || res_name == "G" || res_name == "U" || res_name == "C" ||
                res_name == "DA" || res_name == "DG" || res_name == "DT" || res_name == "DC") {
                continue;
            }
            if (is_het_residue_gemmi(residue_p)) continue;
            residue_vec.push_back(residue_p);
        }
        if (static_cast<int>(residue_vec.size()) >= static_cast<int>(n_residues_per_blob)) {
            if (!residue_vec.empty()) {
                residue_run_t_gemmi rr(n_residues_per_blob);
                rr.add(residue_vec);
                residue_runs.push_back(rr);
                residue_vec.clear();
            }
        }
    }

    // Calculate atom map from the model
    clipper::Xmap<float> calc_map = calc_atom_map_from_model(
        model, xmap.cell(), xmap.spacegroup(), xmap.grid_sampling());

    if (calc_map.is_null()) {
        std::cout << "OOPS! calc_map is null" << std::endl;
        return std::make_pair(res_map_all_atom, res_map_side_chain);
    }

    // Build contributor map - map each grid point to the set of residues that contribute
    clipper::Xmap<std::set<const gemmi::Residue *>> contributor_map(
        xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());

    float atom_radius_sq = atom_mask_radius * atom_mask_radius;
    for (const gemmi::Chain &chain : model.chains) {
        for (const gemmi::Residue &residue : chain.residues) {
            for (const gemmi::Atom &at : residue.atoms) {
                clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
                clipper::Coord_frac cf = co.coord_frac(contributor_map.cell());
                clipper::Coord_frac box0(
                    cf.u() - atom_mask_radius / contributor_map.cell().descr().a(),
                    cf.v() - atom_mask_radius / contributor_map.cell().descr().b(),
                    cf.w() - atom_mask_radius / contributor_map.cell().descr().c());
                clipper::Coord_frac box1(
                    cf.u() + atom_mask_radius / contributor_map.cell().descr().a(),
                    cf.v() + atom_mask_radius / contributor_map.cell().descr().b(),
                    cf.w() + atom_mask_radius / contributor_map.cell().descr().c());
                clipper::Grid_map grid(
                    box0.coord_grid(contributor_map.grid_sampling()),
                    box1.coord_grid(contributor_map.grid_sampling()));
                clipper::Xmap_base::Map_reference_coord ix(contributor_map, grid.min()), iu, iv, iw;
                for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u()) {
                    for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v()) {
                        for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w()) {
                            if ((iw.coord().coord_frac(contributor_map.grid_sampling())
                                    .coord_orth(contributor_map.cell()) - co).lengthsq() < atom_radius_sq) {
                                contributor_map[iw].insert(&residue);
                            }
                        }
                    }
                }
            }
        }
    }

    // Get stats for each residue run
    for (const residue_run_t_gemmi &residue_run : residue_runs) {
        coot::util::density_correlation_stats_info_t dcs_all_atom;
        coot::util::density_correlation_stats_info_t dcs_side_chain;
        std::vector<clipper::Coord_orth> atom_positions_to_avoid;

        if (exclude_CON) {
            for (const gemmi::Residue *residue_p : residue_run.residues) {
                for (const gemmi::Atom &at : residue_p->atoms) {
                    // Use trimmed comparison - CIF files don't pad atom names
                    if (atom_name_is(at, "N")) atom_positions_to_avoid.push_back(coot::co_gemmi(&at));
                    if (atom_name_is(at, "C")) atom_positions_to_avoid.push_back(coot::co_gemmi(&at));
                    if (atom_name_is(at, "O")) atom_positions_to_avoid.push_back(coot::co_gemmi(&at));
                    if (atom_name_is(at, "H")) atom_positions_to_avoid.push_back(coot::co_gemmi(&at));
                }
            }
        }

        float atom_radius = atom_mask_radius;
        for (const gemmi::Residue *residue_p : residue_run.residues) {
            for (const gemmi::Atom &at : residue_p->atoms) {
                clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
                clipper::Coord_frac cf = co.coord_frac(contributor_map.cell());
                clipper::Coord_frac box0(
                    cf.u() - atom_radius / contributor_map.cell().descr().a(),
                    cf.v() - atom_radius / contributor_map.cell().descr().b(),
                    cf.w() - atom_radius / contributor_map.cell().descr().c());
                clipper::Coord_frac box1(
                    cf.u() + atom_radius / contributor_map.cell().descr().a(),
                    cf.v() + atom_radius / contributor_map.cell().descr().b(),
                    cf.w() + atom_radius / contributor_map.cell().descr().c());
                clipper::Grid_map grid(
                    box0.coord_grid(contributor_map.grid_sampling()),
                    box1.coord_grid(contributor_map.grid_sampling()));
                float a_rad_sq = atom_radius * atom_radius;
                clipper::Xmap_base::Map_reference_coord ix(contributor_map, grid.min()), iu, iv, iw;
                for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u()) {
                    for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v()) {
                        for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w()) {
                            if ((iw.coord().coord_frac(contributor_map.grid_sampling())
                                    .coord_orth(contributor_map.cell()) - co).lengthsq() < a_rad_sq) {
                                bool found_other_residue_flag = false;
                                for (const gemmi::Residue *contributor_residue : contributor_map[iw]) {
                                    bool in_run = false;
                                    for (const gemmi::Residue *rr : residue_run.residues) {
                                        if (rr == contributor_residue) {
                                            in_run = true;
                                            break;
                                        }
                                    }
                                    if (!in_run) {
                                        found_other_residue_flag = true;
                                        break;
                                    }
                                }
                                if (!found_other_residue_flag) {
                                    const float &xf = xmap[iw];
                                    const float &yf = calc_map[iw];
                                    if (!clipper::Util::is_nan(yf)) {
                                        double x = static_cast<double>(xf);
                                        double y = static_cast<double>(yf);
                                        if (exclude_CON) {
                                            bool found_a_close_NOC = false;
                                            const double dd_inside = NOC_mask_radius * NOC_mask_radius;
                                            for (const auto &NOC_pos : atom_positions_to_avoid) {
                                                double d = (co - NOC_pos).lengthsq();
                                                if (d < dd_inside) {
                                                    found_a_close_NOC = true;
                                                    break;
                                                }
                                            }
                                            if (found_a_close_NOC)
                                                dcs_all_atom.add(x, y);
                                            else
                                                dcs_side_chain.add(x, y);
                                        } else {
                                            dcs_all_atom.add(x, y);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        const gemmi::Residue *mid_res = residue_run.residue_mid();
        if (mid_res) {
            coot::residue_spec_t spec(chain_id,
                                      mid_res->seqid.num.value,
                                      std::string(1, mid_res->seqid.icode));
            res_map_all_atom[spec] = dcs_all_atom;
            res_map_side_chain[spec] = dcs_side_chain;
        }
    }

    return std::make_pair(res_map_all_atom, res_map_side_chain);
}

} // namespace util
} // namespace coot