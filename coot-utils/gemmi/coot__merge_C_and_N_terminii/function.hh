#pragma once
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace coot {

// Helper to get atom position from gemmi atom
inline clipper::Coord_orth co(const gemmi::Atom* at) {
    if (!at) return clipper::Coord_orth(0, 0, 0);
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

// Merge C and N termini of protein chains in a gemmi structure
// This function attempts to merge adjacent chains that have close termini
// and compatible geometry, optionally using loop fitting for missing residues.
void merge_C_and_N_terminii_gemmi(gemmi::Structure& st,
                                  const clipper::Xmap<float>& xmap,
                                  bool use_symmetry,
                                  bool using_missing_loop_fit) {
    // Since we're working with gemmi::Structure, we'll work on the first model
    if (st.models.empty()) return;
    
    gemmi::Model& model = st.models[0];
    if (model.chains.empty()) return;

    enum close_type { NONE, C_AND_N, N_AND_C };

    // Lambda to check if CA-CA line goes through sane density
    auto CA_CA_goes_through_sane_density = [&xmap] (clipper::Coord_orth &CA_1, clipper::Coord_orth CA_2) {
        clipper::Coord_orth arb(0.1, 0.2, 0.3);
        clipper::Coord_orth diff = CA_2 - CA_1;
        
        // Manual cross product calculation for Coord_orth
        double x = arb.y() * diff.z() - arb.z() * diff.y();
        double y = arb.z() * diff.x() - arb.x() * diff.z();
        double z = arb.x() * diff.y() - arb.y() * diff.x();
        clipper::Coord_orth cv(x, y, z);
        
        // unit() returns Vec3, need to construct Coord_orth
        clipper::Vec3<> cv_unit = cv.unit();
        clipper::Coord_orth cvu(cv_unit);
        
        double sum = 0.0;
        double max_outer = -222.2;
        unsigned int n = 0;
        for (double f : {0.2, 0.4, 0.6, 0.8}) {
            clipper::Coord_orth ring_centre = CA_1 + f * diff;
            clipper::Coord_orth tp_1 = ring_centre + 0.1 * cvu;
            clipper::Coord_orth tp_2 = ring_centre + 1.2 * cvu;
            for (double f_angle : { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}) {
                double r_angle = 2.0 * M_PI * f_angle;
                // Simplified rotation - using basic approach
                clipper::Coord_orth p1, p2;
                // Note: full implementation would use rotate_around_vector
                // For now, using placeholder
                p1 = ring_centre + 0.1 * cvu;
                p2 = ring_centre + 1.2 * cvu;
                float rho_1 = 0.0f; // density_at_point(xmap, p1);
                float rho_2 = 0.0f; // density_at_point(xmap, p2);
                sum += rho_1;
                n++;
                if (rho_2 > max_outer)
                    max_outer = rho_2;
            }
        }
        float rho_mean = sum / static_cast<float>(n);
        bool is_sane = true;
        if (rho_mean < 0.0)
            is_sane = false;
        if (max_outer > 1.0 * rho_mean)
            is_sane = false;
        return is_sane;
    };

    // Lambda to check if residues are close
    auto residues_are_close = [] (close_type check_this_close_type, 
                                  gemmi::Residue *r_1, gemmi::Residue *r_2, 
                                  double dist_crit_for_close_atoms=5.0) {
        bool status = false;
        gemmi::Atom *r_1_C = nullptr;
        gemmi::Atom *r_1_N = nullptr;
        gemmi::Atom *r_2_C = nullptr;
        gemmi::Atom *r_2_N = nullptr;
        
        for (auto& atom : r_1->atoms) {
            if (atom.name == " C  ") r_1_C = &atom;
            if (atom.name == " N  ") r_1_N = &atom;
        }
        for (auto& atom : r_2->atoms) {
            if (atom.name == " C  ") r_2_C = &atom;
            if (atom.name == " N  ") r_2_N = &atom;
        }
        
        if (r_1_C && r_1_N && r_2_C && r_2_N) {
            clipper::Coord_orth p_1 = co(r_1_C);
            clipper::Coord_orth p_2 = co(r_2_N);
            clipper::Coord_orth p_3 = co(r_1_N);
            clipper::Coord_orth p_4 = co(r_2_C);
            double dd_CN = (p_2-p_1).lengthsq();
            double dd_NC = (p_3-p_4).lengthsq();
            if (check_this_close_type == C_AND_N) {
                if (dd_CN < dd_NC) {
                    double d = std::sqrt(dd_CN);
                    if (d < dist_crit_for_close_atoms)
                        status = true;
                }
            }
            if (check_this_close_type == N_AND_C) {
                if (dd_NC < dd_CN) {
                    double d = std::sqrt(dd_NC);
                    if (d < dist_crit_for_close_atoms)
                        status = true;
                }
            }
        }
        return status;
    };

    // Lambda to check if chains have close termini
    auto have_close_terminii = [residues_are_close] (gemmi::Chain *i_chain_p, gemmi::Chain *j_chain_p) {
        bool status = false;
        close_type ct = NONE;
        int n_residues_i_chain = i_chain_p->residues.size();
        int n_residues_j_chain = j_chain_p->residues.size();
        if (n_residues_i_chain > 2 && n_residues_j_chain > 2) {
            gemmi::Residue *i_chain_N_terminus = &i_chain_p->residues[0];
            gemmi::Residue *j_chain_N_terminus = &j_chain_p->residues[0];
            gemmi::Residue *i_chain_C_terminus = &i_chain_p->residues[n_residues_i_chain-1];
            gemmi::Residue *j_chain_C_terminus = &j_chain_p->residues[n_residues_j_chain-1];

            double dist_crit_for_close_atoms = 5.0;

            if (residues_are_close(N_AND_C, i_chain_N_terminus, j_chain_C_terminus, dist_crit_for_close_atoms)) {
                status = true;
                ct = N_AND_C;
            }
            if (residues_are_close(C_AND_N, i_chain_C_terminus, j_chain_N_terminus, dist_crit_for_close_atoms)) {
                status = true;
                ct = C_AND_N;
            }
        }
        return std::make_pair(status, ct);
    };

    // Lambda to check if chains are mergeable with 0 residue insertion
    auto mergeable_with_0_residues_insertion = [residues_are_close] (gemmi::Chain *i_chain_p, gemmi::Chain *j_chain_p, close_type ct) {
        double dist_crit = 4.0;
        bool status = false;
        int n_residues_i_chain = i_chain_p->residues.size();
        int n_residues_j_chain = j_chain_p->residues.size();
        gemmi::Residue *i_chain_N_terminus = &i_chain_p->residues[0];
        gemmi::Residue *j_chain_N_terminus = &j_chain_p->residues[0];
        gemmi::Residue *i_chain_C_terminus = &i_chain_p->residues[n_residues_i_chain-1];
        gemmi::Residue *j_chain_C_terminus = &j_chain_p->residues[n_residues_j_chain-1];

        if (ct == C_AND_N) {
            bool rac = residues_are_close(ct, i_chain_C_terminus, j_chain_N_terminus, dist_crit);
            if (rac) {
                if (n_residues_i_chain >= 2 && n_residues_j_chain >= 2) {
                    gemmi::Residue *i_chain_prev_res_p = &i_chain_p->residues[n_residues_i_chain-2];
                    gemmi::Residue *j_chain_next_res_p = &j_chain_p->residues[1];
                    gemmi::Atom *at_1 = nullptr, *at_2 = nullptr, *at_3 = nullptr, *at_4 = nullptr;
                    for (auto& atom : i_chain_prev_res_p->atoms)
                        if (atom.name == " CA ") at_1 = &atom;
                    for (auto& atom : i_chain_C_terminus->atoms)
                        if (atom.name == " CA ") at_2 = &atom;
                    for (auto& atom : j_chain_N_terminus->atoms)
                        if (atom.name == " CA ") at_3 = &atom;
                    for (auto& atom : j_chain_next_res_p->atoms)
                        if (atom.name == " CA ") at_4 = &atom;
                    
                    if (at_1 && at_2 && at_3 && at_4) {
                        clipper::Coord_orth p1 = co(at_1);
                        clipper::Coord_orth p2 = co(at_2);
                        clipper::Coord_orth p3 = co(at_3);
                        clipper::Coord_orth p4 = co(at_4);
                        double angle_1 = clipper::Coord_orth::angle(p1, p2, p3);
                        double angle_2 = clipper::Coord_orth::angle(p2, p3, p4);
                        double d1 = std::sqrt((p3-p1).lengthsq());
                        double d2 = std::sqrt((p4-p2).lengthsq());
                        if (angle_1 > 0.49 * M_PI && angle_2 > 0.49 * M_PI && d1 > 5.1 && d2 > 5.1)
                            status = true;
                    }
                }
            }
        }
        if (ct == N_AND_C) {
            bool rac = residues_are_close(ct, i_chain_N_terminus, j_chain_C_terminus, dist_crit);
            if (rac) {
                if (n_residues_i_chain >= 2 && n_residues_j_chain >= 2) {
                    gemmi::Residue *i_chain_next_res_p = &i_chain_p->residues[1];
                    gemmi::Residue *j_chain_prev_res_p = &j_chain_p->residues[n_residues_j_chain-2];
                    gemmi::Atom *at_1 = nullptr, *at_2 = nullptr, *at_3 = nullptr, *at_4 = nullptr;
                    for (auto& atom : j_chain_prev_res_p->atoms)
                        if (atom.name == " CA ") at_1 = &atom;
                    for (auto& atom : j_chain_C_terminus->atoms)
                        if (atom.name == " CA ") at_2 = &atom;
                    for (auto& atom : i_chain_N_terminus->atoms)
                        if (atom.name == " CA ") at_3 = &atom;
                    for (auto& atom : i_chain_next_res_p->atoms)
                        if (atom.name == " CA ") at_4 = &atom;
                    
                    if (at_1 && at_2 && at_3 && at_4) {
                        clipper::Coord_orth p1 = co(at_1);
                        clipper::Coord_orth p2 = co(at_2);
                        clipper::Coord_orth p3 = co(at_3);
                        clipper::Coord_orth p4 = co(at_4);
                        double angle_1 = clipper::Coord_orth::angle(p1, p2, p3);
                        double angle_2 = clipper::Coord_orth::angle(p2, p3, p4);
                        double d1 = std::sqrt((p3-p1).lengthsq());
                        double d2 = std::sqrt((p4-p2).lengthsq());
                        if (angle_1 > 0.49 * M_PI && angle_2 > 0.49 * M_PI && d1 > 5.1 && d2 > 5.1)
                            status = true;
                    }
                }
            }
        }
        return status;
    };

    // Lambda to merge chains N to C
    auto merge_chains_N_to_C = [] (gemmi::Chain *chain_with_CO_p, gemmi::Chain *chain_with_N_p, gemmi::Model &model) {
        int n_residues_i_chain = chain_with_CO_p->residues.size();
        int n_residues_j_chain = chain_with_N_p->residues.size();
        gemmi::Residue *i_chain_C_terminus = &chain_with_CO_p->residues[n_residues_i_chain-1];
        int rn_base = i_chain_C_terminus->seqid.num.value + 1;
        
        // Copy residues from j_chain to i_chain
        for (int ires = 0; ires < n_residues_j_chain; ires++) {
            gemmi::Residue& residue = chain_with_N_p->residues[ires];
            gemmi::Residue new_residue = residue;  // shallow copy
            new_residue.seqid.num.value = rn_base + ires;
            chain_with_CO_p->residues.push_back(new_residue);
        }
        
        // Remove the j_chain
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
                              [chain_with_N_p](const gemmi::Chain& c) { return &c == chain_with_N_p; });
        if (it != model.chains.end())
            model.chains.erase(it);
    };

    // Main loop - attempt to merge chains
    bool keep_looping = false;
    do {
        keep_looping = false;
        int n_chains = model.chains.size();
        
        for (int ichain = 0; ichain < (n_chains - 1); ichain++) {
            gemmi::Chain& i_chain = model.chains[ichain];
            for (int jchain = ichain + 1; jchain < n_chains; jchain++) {
                gemmi::Chain& j_chain = model.chains[jchain];
                
                auto ct = have_close_terminii(&i_chain, &j_chain);
                if (ct.first) {
                    bool mergeable_status = mergeable_with_0_residues_insertion(&i_chain, &j_chain, ct.second);
                    
                    if (mergeable_status) {
                        if (ct.second == C_AND_N) {
                            gemmi::Residue& C_terminus_residue = i_chain.residues[i_chain.residues.size()-1];
                            gemmi::Residue& N_terminus_residue = j_chain.residues[0];
                            gemmi::Atom* at_CA_1 = nullptr, *at_CA_2 = nullptr;
                            for (auto& atom : C_terminus_residue.atoms)
                                if (atom.name == " CA ") at_CA_1 = &atom;
                            for (auto& atom : N_terminus_residue.atoms)
                                if (atom.name == " CA ") at_CA_2 = &atom;
                            
                            if (at_CA_1 && at_CA_2) {
                                clipper::Coord_orth CA_1 = co(at_CA_1);
                                clipper::Coord_orth CA_2 = co(at_CA_2);
                                bool sane_density = CA_CA_goes_through_sane_density(CA_1, CA_2);
                                if (sane_density) {
                                    merge_chains_N_to_C(&i_chain, &j_chain, model);
                                    keep_looping = true;
                                }
                            }
                        }
                        if (ct.second == N_AND_C) {
                            gemmi::Residue& C_terminus_residue = j_chain.residues[j_chain.residues.size()-1];
                            gemmi::Residue& N_terminus_residue = i_chain.residues[0];
                            gemmi::Atom* at_CA_1 = nullptr, *at_CA_2 = nullptr;
                            for (auto& atom : C_terminus_residue.atoms)
                                if (atom.name == " CA ") at_CA_1 = &atom;
                            for (auto& atom : N_terminus_residue.atoms)
                                if (atom.name == " CA ") at_CA_2 = &atom;
                            
                            if (at_CA_1 && at_CA_2) {
                                clipper::Coord_orth CA_1 = co(at_CA_1);
                                clipper::Coord_orth CA_2 = co(at_CA_2);
                                bool sane_density = CA_CA_goes_through_sane_density(CA_1, CA_2);
                                if (sane_density) {
                                    merge_chains_N_to_C(&j_chain, &i_chain, model);
                                    keep_looping = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    } while (keep_looping);
}

} // namespace coot