#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <coot/coords/Cartesian.hh>
#include <coot/coords/loop-path.hh>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__loop_path/gemmi/function.hh"

namespace coot {

/**
 * gemmi port of Bond_lines_container::do_Ca_loop
 *
 * Represents missing residues as a curved loop between four anchor atoms.
 * Returns the loop path segments and a "bad gap" flag.
 */
std::pair<bool, std::vector<coot::CartesianPair>> do_Ca_loop_gemmi(
    int ires,
    int nres,
    const gemmi::Chain* chain_p,
    const gemmi::Residue* residue_prev,
    const gemmi::Residue* residue_this) {

    std::pair<bool, std::vector<coot::CartesianPair>> result{false, {}};

    if (!chain_p || !residue_prev || !residue_this) return result;

    std::string res_name_1 = residue_prev->name;
    std::string res_name_2 = residue_this->name;

    if (res_name_1 == "HOH") return result;
    if (res_name_2 == "HOH") return result;

    int res_idx_n_start_back = (ires - 1) - 2;
    int res_idx_n_end_forwards = ires + 2;
    if (res_idx_n_start_back < 0) return result;
    if (res_idx_n_end_forwards >= nres) return result;

    // Find the residues at these indices
    const gemmi::Residue* res_start_back = nullptr;
    const gemmi::Residue* res_end_forwards = nullptr;
    {
        unsigned int idx = 0;
        for (const auto& res : chain_p->residues) {
            if (static_cast<int>(idx) == res_idx_n_start_back) res_start_back = &res;
            if (static_cast<int>(idx) == res_idx_n_end_forwards) res_end_forwards = &res;
            idx++;
        }
    }

    if (!res_start_back || !res_end_forwards) return result;

    // Helper: check if atom name matches (handles both PDB-padded and CIF-unpadded)
    auto matches_atom = [](const std::string& an, const std::string& target) {
        if (an == target) return true;
        if (target.size() <= 4) {
            std::string padded = target;
            padded.resize(4, ' ');
            if (an == padded) return true;
        }
        return false;
    };

    // Find CA or P atoms in each of the 4 residues
    const gemmi::Atom* at_pp_1 = nullptr; // res_start_back (n-2)
    const gemmi::Atom* at_pp_2 = nullptr; // residue_prev (n-1)
    const gemmi::Atom* at_pp_3 = nullptr; // residue_this (n)
    const gemmi::Atom* at_pp_4 = nullptr; // res_end_forwards (n+2)

    // Find CA or P in res_start_back
    for (const auto& a : res_start_back->atoms) {
        if (matches_atom(a.name, "CA") || matches_atom(a.name, "P")) {
            at_pp_1 = &a;
            break;
        }
    }
    if (!at_pp_1) return result;

    // Find CA, P, or O3' in residue_prev (prefer O3')
    for (const auto& a : residue_prev->atoms) {
        if (matches_atom(a.name, "CA") || matches_atom(a.name, "P")) {
            at_pp_2 = &a;
        }
        if (matches_atom(a.name, "O3'")) {
            at_pp_2 = &a;
            break; // Prefer O3'
        }
    }
    if (!at_pp_2) return result;

    // Find CA or P in residue_this
    for (const auto& a : residue_this->atoms) {
        if (matches_atom(a.name, "CA") || matches_atom(a.name, "P")) {
            at_pp_3 = &a;
            break;
        }
    }
    if (!at_pp_3) return result;

    // Find CA or P in res_end_forwards
    for (const auto& a : res_end_forwards->atoms) {
        if (matches_atom(a.name, "CA") || matches_atom(a.name, "P")) {
            at_pp_4 = &a;
            break;
        }
    }
    if (!at_pp_4) return result;

    // Are C of previous and N of next close as in a peptide bond?
    {
        const gemmi::Atom* C_prev = nullptr;
        for (const auto& a : residue_prev->atoms) {
            if (matches_atom(a.name, "C")) { C_prev = &a; break; }
        }
        if (C_prev) {
            const gemmi::Atom* N_this = nullptr;
            for (const auto& a : residue_this->atoms) {
                if (matches_atom(a.name, "N")) { N_this = &a; break; }
            }
            if (N_this) {
                double dx = C_prev->pos.x - N_this->pos.x;
                double dy = C_prev->pos.y - N_this->pos.y;
                double dz = C_prev->pos.z - N_this->pos.z;
                float dist_sqrd = static_cast<float>(dx * dx + dy * dy + dz * dz);
                if (dist_sqrd < 2.5f * 2.5f)
                    return result; // C and N are close — no loop
            }
        }
    }

    // Are O3' of previous and P of next close as in a phosphodiester?
    {
        const gemmi::Atom* O3prime_prev = nullptr;
        for (const auto& a : residue_prev->atoms) {
            if (matches_atom(a.name, "O3'")) { O3prime_prev = &a; break; }
        }
        if (O3prime_prev) {
            const gemmi::Atom* P_this = nullptr;
            for (const auto& a : residue_this->atoms) {
                if (matches_atom(a.name, "P")) { P_this = &a; break; }
            }
            if (P_this) {
                double dx = O3prime_prev->pos.x - P_this->pos.x;
                double dy = O3prime_prev->pos.y - P_this->pos.y;
                double dz = O3prime_prev->pos.z - P_this->pos.z;
                float dist_sqrd = static_cast<float>(dx * dx + dy * dy + dz * dz);
                if (dist_sqrd < 2.5f * 2.5f)
                    return result; // O3' and P are close — no loop
            }
        }
    }

    // Build CRAs for the four atoms
    gemmi::CRA cra_back_2 = { const_cast<gemmi::Chain*>(chain_p),
                              const_cast<gemmi::Residue*>(res_start_back),
                              const_cast<gemmi::Atom*>(at_pp_1) };
    gemmi::CRA cra_prev = { const_cast<gemmi::Chain*>(chain_p),
                            const_cast<gemmi::Residue*>(residue_prev),
                            const_cast<gemmi::Atom*>(at_pp_2) };
    gemmi::CRA cra_this = { const_cast<gemmi::Chain*>(chain_p),
                            const_cast<gemmi::Residue*>(residue_this),
                            const_cast<gemmi::Atom*>(at_pp_3) };
    gemmi::CRA cra_end = { const_cast<gemmi::Chain*>(chain_p),
                           const_cast<gemmi::Residue*>(res_end_forwards),
                           const_cast<gemmi::Atom*>(at_pp_4) };

    coot::Cartesian pp_2(at_pp_2->pos.x, at_pp_2->pos.y, at_pp_2->pos.z);
    coot::Cartesian pp_3(at_pp_3->pos.x, at_pp_3->pos.y, at_pp_3->pos.z);
    float a = (pp_3 - pp_2).amplitude();
    int n_line_segments = static_cast<int>(a * 1.2);

    return coot::loop_path_gemmi(cra_back_2, cra_prev, cra_this, cra_end,
                                 static_cast<unsigned int>(n_line_segments));
}

} // namespace coot
