#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__cablam_markup_t__cablam_markup_t/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

inline cablam_markup_t_gemmi calc_cablam_gemmi(
    gemmi::Chain *chain_p,
    gemmi::Residue *residue_this_p,
    int ires,
    double score)
{
    cablam_markup_t_gemmi cm; // default: score=-1, residue=nullptr

    if (ires - 1 < 0 || ires + 1 >= static_cast<int>(chain_p->residues.size())) {
        return cm;
    }

    gemmi::Residue *residue_prev_p = &chain_p->residues[ires - 1];
    gemmi::Residue *residue_next_p = &chain_p->residues[ires + 1];

    // tandem residues checks
    if (residue_prev_p->seqid.num.value + 1 != residue_this_p->seqid.num.value) {
        return cm;
    }
    if (residue_next_p->seqid.num.value - 1 != residue_this_p->seqid.num.value) {
        return cm;
    }

    // Find relevant atoms
    gemmi::Atom *O_prev = nullptr;
    gemmi::Atom *O_this = nullptr;
    gemmi::Atom *O_next = nullptr;
    gemmi::Atom *CA_prev = nullptr;
    gemmi::Atom *CA_this = nullptr;
    gemmi::Atom *CA_next = nullptr;
    gemmi::Atom *CA_next_next = nullptr;

    for (gemmi::Atom &at : residue_prev_p->atoms) {
        if (at.altloc == '\0') {
            if (at.name == "O")  { O_prev = &at; }
            else if (at.name == "CA") { CA_prev = &at; }
        }
    }
    for (gemmi::Atom &at : residue_this_p->atoms) {
        if (at.altloc == '\0') {
            if (at.name == "O")  { O_this = &at; }
            else if (at.name == "CA") { CA_this = &at; }
        }
    }
    for (gemmi::Atom &at : residue_next_p->atoms) {
        if (at.altloc == '\0') {
            if (at.name == "O")  { O_next = &at; }
            else if (at.name == "CA") { CA_next = &at; }
        }
    }
    if (ires + 2 < static_cast<int>(chain_p->residues.size())) {
        for (gemmi::Atom &at : chain_p->residues[ires + 2].atoms) {
            if (at.altloc == '\0') {
                if (at.name == "CA") { CA_next_next = &at; }
            }
        }
    }

    if (!O_prev || !O_this || !O_next) return cm;
    if (!CA_prev || !CA_this || !CA_next || !CA_next_next) return cm;

    // Check CA-CA distances
    clipper::Coord_orth p1 = coot::co_gemmi(CA_prev);
    clipper::Coord_orth p2 = coot::co_gemmi(CA_this);
    clipper::Coord_orth p3 = coot::co_gemmi(CA_next);

    clipper::Coord_orth diff1 = p2 - p1;
    clipper::Coord_orth diff2 = p2 - p3;
    double v1_sqrd = diff1[0]*diff1[0] + diff1[1]*diff1[1] + diff1[2]*diff1[2];
    double v2_sqrd = diff2[0]*diff2[0] + diff2[1]*diff2[1] + diff2[2]*diff2[2];
    double v1 = std::sqrt(v1_sqrd);
    double v2 = std::sqrt(v2_sqrd);

    if (v1 >= 3.9 || v2 >= 3.9) return cm;

    // Build result — replicate the 7-atom constructor logic from cablam-markup.cc
    cablam_markup_t_gemmi cm_local;
    cm_local.residue = residue_this_p;
    cm_local.score = -1;

    cm_local.O_prev_pos = coot::co_gemmi(O_prev);
    cm_local.O_this_pos = coot::co_gemmi(O_this);
    cm_local.O_next_pos = coot::co_gemmi(O_next);

    // CA projection point: "this" pair (CA_this -> CA_next)
    {
        clipper::Coord_orth d = cm_local.O_this_pos - p2;
        clipper::Coord_orth r = p3 - p2;
        double r_len_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        double r_len = std::sqrt(r_len_sq);
        if (r_len > 1e-10) {
            double dp = r[0]*d[0] + r[1]*d[1] + r[2]*d[2];
            double r1_len = dp / r_len;
            cm_local.CA_proj_point_this = p2 + r1_len * clipper::Coord_orth(r[0]/r_len, r[1]/r_len, r[2]/r_len);
        } else {
            cm_local.CA_proj_point_this = p2;
        }
    }

    // CA projection point: "prev" pair (CA_prev -> CA_this)
    {
        clipper::Coord_orth d = cm_local.O_prev_pos - p1;
        clipper::Coord_orth r = p2 - p1;
        double r_len_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        double r_len = std::sqrt(r_len_sq);
        if (r_len > 1e-10) {
            double dp = r[0]*d[0] + r[1]*d[1] + r[2]*d[2];
            double r1_len = dp / r_len;
            cm_local.CA_proj_point_prev = p1 + r1_len * clipper::Coord_orth(r[0]/r_len, r[1]/r_len, r[2]/r_len);
        } else {
            cm_local.CA_proj_point_prev = p1;
        }
    }

    // CA projection point: "next" pair (CA_next -> CA_next_next)
    {
        clipper::Coord_orth CA_next_next_pos = coot::co_gemmi(CA_next_next);
        clipper::Coord_orth d = cm_local.O_next_pos - p3;
        clipper::Coord_orth r = CA_next_next_pos - p3;
        double r_len_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        double r_len = std::sqrt(r_len_sq);
        if (r_len > 1e-10) {
            double dp = r[0]*d[0] + r[1]*d[1] + r[2]*d[2];
            double r1_len = dp / r_len;
            cm_local.CA_proj_point_next = p3 + r1_len * clipper::Coord_orth(r[0]/r_len, r[1]/r_len, r[2]/r_len);
        } else {
            cm_local.CA_proj_point_next = p3;
        }
    }

    cm_local.score = score;
    cm = cm_local;

    return cm;
}

} // namespace coot