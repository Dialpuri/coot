#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <vector>
#include <string>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__lidia_utils__get_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__lidia_utils__co/gemmi/function.hh"

namespace coot {
namespace chem_feat_clust {

// Gemmi-adapted water attribute: replaces mmdb::Atom* with gemmi::Atom*
struct water_attribs_gemmi {
    int imol;
    unsigned int ligand_idx;
    unsigned int water_spec_idx;
    gemmi::Atom* atom_p;
    clipper::Coord_orth co;

    water_attribs_gemmi(int i, unsigned l, unsigned w, gemmi::Atom* a, clipper::Coord_orth c)
        : imol(i), ligand_idx(l), water_spec_idx(w), atom_p(a), co(c) {}
};

// Gemmi-adapted solvated ligand spec: replaces mmdb::Manager* with gemmi::Structure*
struct chem_feat_solvated_ligand_spec_gemmi {
    coot::residue_spec_t ligand_spec;
    std::vector<coot::residue_spec_t> waters;
    gemmi::Structure* mol;
    int imol;
};

// Gemmi-adapted cluster
struct chem_feat_clust_gemmi {
    std::vector<coot::residue_spec_t> protein_residues;
    std::vector<chem_feat_solvated_ligand_spec_gemmi> ligands;
    std::vector<water_attribs_gemmi> water_positions;
    bool setup_success;
    double water_dist_cutoff;

    chem_feat_clust_gemmi(
        const std::vector<coot::residue_spec_t>& prot_res,
        const std::vector<chem_feat_solvated_ligand_spec_gemmi>& lig_spec,
        double wd_cut,
        bool success = true)
        : protein_residues(prot_res),
          ligands(lig_spec),
          water_dist_cutoff(wd_cut),
          setup_success(success) {}
};

// Check whether pt is within dist_cutoff of any coordinate in all_ligand_coords
inline bool is_close_to_a_ligand_atom(
    const clipper::Coord_orth& pt,
    const std::vector<clipper::Coord_orth>& all_ligand_coords,
    double dist_cutoff)
{
    for (const auto& lc : all_ligand_coords) {
        clipper::Coord_orth diff(pt - lc);
        double dist = std::sqrt(diff.x() * diff.x() +
                                diff.y() * diff.y() +
                                diff.z() * diff.z());
        if (dist < dist_cutoff) return true;
    }
    return false;
}

// Gemmi port of fill_waters: iterates over ligand waters, finds HOH O atoms
// close to ligand atoms, and records them in water_positions.
inline void fill_waters_gemmi(chem_feat_clust_gemmi& clust)
{
    // Build a vector of all ligand atom coordinates
    std::vector<clipper::Coord_orth> lig_coords;
    for (const auto& ligand : clust.ligands) {
        gemmi::Residue* res = coot::lidia_utils::get_residue_gemmi(
            ligand.ligand_spec, ligand.mol);
        if (res) {
            for (const gemmi::Atom& atom : res->atoms) {
                lig_coords.push_back(coot::lidia_utils::co_gemmi(&atom));
            }
        }
    }

    // Iterate over every water spec attached to each ligand
    for (unsigned int ilig = 0; ilig < clust.ligands.size(); ++ilig) {
        for (unsigned int iw = 0; iw < clust.ligands[ilig].waters.size(); ++iw) {
            gemmi::Residue* res = coot::lidia_utils::get_residue_gemmi(
                clust.ligands[ilig].waters[iw], clust.ligands[ilig].mol);
            if (!res) continue;

            if (res->name == "HOH") {
                // Find the O atom inside this HOH residue
                const gemmi::Atom* o_atom = nullptr;
                for (const gemmi::Atom& atom : res->atoms) {
                    // Trim whitespace from the atom name and compare
                    std::string trimmed = atom.name;
                    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
                    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
                    if (trimmed == "O") {
                        o_atom = &atom;
                        break;
                    }
                }

                if (o_atom) {
                    clipper::Coord_orth pt = coot::lidia_utils::co_gemmi(o_atom);

                    if (is_close_to_a_ligand_atom(pt, lig_coords, clust.water_dist_cutoff)) {
                        water_attribs_gemmi wa(
                            clust.ligands[ilig].imol,
                            ilig,
                            iw,
                            const_cast<gemmi::Atom*>(o_atom),
                            pt);
                        clust.water_positions.push_back(wa);
                    }
                }
            }
        }
    }
}

} // namespace chem_feat_clust
} // namespace coot