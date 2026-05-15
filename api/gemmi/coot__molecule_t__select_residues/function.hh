#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <optional>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

// Include ports for coot utilities we call
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_following_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_previous_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__residues_in_chain/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residues_near_residue/gemmi/function.hh"

namespace coot {

// Helper: convert a gemmi::Residue* to residue_spec_t by searching in the structure.
static inline residue_spec_t
spec_from_residue(gemmi::Residue *res, const gemmi::Structure &st) {
    if (!res) return residue_spec_t();
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (&r == res) {
                    std::string ic = std::string(1, r.seqid.icode);
                    if (ic == " ") ic = "";
                    return residue_spec_t(chain.name, r.seqid.num.value, ic);
                }
            }
        }
    }
    return residue_spec_t();
}

static inline std::vector<gemmi::Residue*>
select_residues_gemmi(const residue_spec_t &residue_spec,
                      const std::string &mode_in,
                      gemmi::Structure &st) {

    auto all_residues = [] (gemmi::Structure &structure) {
        std::vector<gemmi::Residue*> rv;
        for (auto& model : structure.models) {
            for (auto& chain : model.chains) {
                for (auto& residue : chain.residues) {
                    rv.push_back(&residue);
                }
            }
        }
        return rv;
    };

    std::string mode = mode_in;
    if (mode == "LITERAL") mode = "SINGLE";

    std::vector<gemmi::Residue*> rv;

    gemmi::Residue *residue_p = util::get_residue_gemmi(
            residue_spec.chain_id,
            residue_spec.res_no,
            residue_spec.ins_code,
            st);

    if (residue_p) {
        if (mode == "SINGLE") {
            rv.push_back(residue_p);
        }
        if (mode == "TRIPLE") {
            std::optional<gemmi::CRA> r_p_1_opt = util::get_following_residue_gemmi(residue_spec, st);

            gemmi::Residue *r_p_1 = r_p_1_opt.has_value() ? r_p_1_opt->residue : nullptr;
            gemmi::Residue *r_m_1 = util::get_previous_residue_gemmi(residue_spec, st);

            if (r_m_1) rv.push_back(r_m_1);
            rv.push_back(residue_p);
            if (r_p_1) rv.push_back(r_p_1);
        }
        if (mode == "QUINTUPLE") {
            std::optional<gemmi::CRA> r_p_1_opt = util::get_following_residue_gemmi(residue_spec, st);
            gemmi::Residue *r_p_1 = r_p_1_opt.has_value() ? r_p_1_opt->residue : nullptr;

            gemmi::Residue *r_m_1 = util::get_previous_residue_gemmi(residue_spec, st);

            gemmi::Residue *r_p_2 = nullptr;
            if (r_p_1) {
                residue_spec_t spec_p1 = spec_from_residue(r_p_1, st);
                std::optional<gemmi::CRA> r_p_2_opt = util::get_following_residue_gemmi(spec_p1, st);
                r_p_2 = r_p_2_opt.has_value() ? r_p_2_opt->residue : nullptr;
            }

            gemmi::Residue *r_m_2 = nullptr;
            if (r_m_1) {
                residue_spec_t spec_m1 = spec_from_residue(r_m_1, st);
                r_m_2 = util::get_previous_residue_gemmi(spec_m1, st);
            }

            if (r_m_2) rv.push_back(r_m_2);
            if (r_m_1) rv.push_back(r_m_1);
            rv.push_back(residue_p);
            if (r_p_1) rv.push_back(r_p_1);
            if (r_p_2) rv.push_back(r_p_2);
        }
        if (mode == "HEPTUPLE") {
            std::optional<gemmi::CRA> r_p_1_opt = util::get_following_residue_gemmi(residue_spec, st);
            gemmi::Residue *r_p_1 = r_p_1_opt.has_value() ? r_p_1_opt->residue : nullptr;

            gemmi::Residue *r_m_1 = util::get_previous_residue_gemmi(residue_spec, st);

            gemmi::Residue *r_p_2 = nullptr;
            if (r_p_1) {
                residue_spec_t spec_p1 = spec_from_residue(r_p_1, st);
                std::optional<gemmi::CRA> r_p_2_opt = util::get_following_residue_gemmi(spec_p1, st);
                r_p_2 = r_p_2_opt.has_value() ? r_p_2_opt->residue : nullptr;
            }

            gemmi::Residue *r_m_2 = nullptr;
            if (r_m_1) {
                residue_spec_t spec_m1 = spec_from_residue(r_m_1, st);
                r_m_2 = util::get_previous_residue_gemmi(spec_m1, st);
            }

            gemmi::Residue *r_p_3 = nullptr;
            if (r_p_2) {
                residue_spec_t spec_p2 = spec_from_residue(r_p_2, st);
                std::optional<gemmi::CRA> r_p_3_opt = util::get_following_residue_gemmi(spec_p2, st);
                r_p_3 = r_p_3_opt.has_value() ? r_p_3_opt->residue : nullptr;
            }

            gemmi::Residue *r_m_3 = nullptr;
            if (r_m_2) {
                residue_spec_t spec_m2 = spec_from_residue(r_m_2, st);
                r_m_3 = util::get_previous_residue_gemmi(spec_m2, st);
            }

            if (r_m_3) rv.push_back(r_m_3);
            if (r_m_2) rv.push_back(r_m_2);
            if (r_m_1) rv.push_back(r_m_1);
            rv.push_back(residue_p);
            if (r_p_1) rv.push_back(r_p_1);
            if (r_p_2) rv.push_back(r_p_2);
            if (r_p_3) rv.push_back(r_p_3);
        }
        if (mode == "CHAIN") {
            rv = util::residues_in_chain_gemmi(st.models[0], residue_spec.chain_id);
        }
        if (mode == "ALL") {
            rv = all_residues(st);
        }
        if (mode == "SPHERE") {
            float radius = 4.2;
            auto v = residues_near_residue_gemmi(residue_spec, st, radius);
            rv.push_back(residue_p);
            for (const auto &rs : v) {
                gemmi::Residue *r = util::get_residue_gemmi(
                        rs.chain_id, rs.res_no, rs.ins_code, st);
                if (r) rv.push_back(r);
            }
        }
        if (mode == "BIG_SPHERE") {
            float radius = 8.0;
            auto v = residues_near_residue_gemmi(residue_spec, st, radius);
            rv.push_back(residue_p);
            for (const auto &rs : v) {
                gemmi::Residue *r = util::get_residue_gemmi(
                        rs.chain_id, rs.res_no, rs.ins_code, st);
                if (r) rv.push_back(r);
            }
        }
    }

    return rv;
}

} // namespace coot