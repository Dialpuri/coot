#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <map>
#include <gemmi/model.hpp>

namespace coot {

class match_container_for_residues_t_gemmi {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::string chain_name_1;
    std::string chain_name_2;
    std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_pairs;

    match_container_for_residues_t_gemmi()
        : residue_1(nullptr), residue_2(nullptr) {}

    match_container_for_residues_t_gemmi(gemmi::Residue *r1, gemmi::Residue *r2,
                                          const std::string& cn1, const std::string& cn2)
        : residue_1(r1), residue_2(r2), chain_name_1(cn1), chain_name_2(cn2) {}

    void add(gemmi::Atom *at_1, gemmi::Atom *at_2) {
        atom_pairs.emplace_back(at_1, at_2);
    }
};

class match_container_t {
public:
    std::vector<match_container_for_residues_t_gemmi> matches;

    void add_gemmi(gemmi::Atom *at_1, gemmi::Atom *at_2, const gemmi::Model &model) {
        gemmi::Residue *res_1 = nullptr;
        gemmi::Residue *res_2 = nullptr;
        std::string cn_1, cn_2;

        for (gemmi::Chain& chain : const_cast<gemmi::Model&>(model).chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (&atom == at_1) { res_1 = &res; cn_1 = chain.name; }
                    if (&atom == at_2) { res_2 = &res; cn_2 = chain.name; }
                }
            }
        }

        if (res_1 && res_2) {
            bool added = false;
            for (size_t i = 0; i < matches.size(); i++) {
                if (matches[i].residue_1 == res_1) {
                    if (matches[i].residue_2 == res_2) {
                        matches[i].add(at_1, at_2);
                        added = true;
                        break;
                    }
                }
            }
            if (!added) {
                match_container_for_residues_t_gemmi m(res_1, res_2, cn_1, cn_2);
                m.add(at_1, at_2);
                matches.push_back(m);
            }
        }
    }

    match_container_for_residues_t_gemmi find_best_match_gemmi() const {
        match_container_for_residues_t_gemmi rm;
        double best_av_devi = 999999999.9;
        for (size_t i = 0; i < matches.size(); i++) {
            const match_container_for_residues_t_gemmi &m = matches[i];
            if (m.atom_pairs.size() > 2) {
                double sum_devi = 0;
                for (size_t iat = 0; iat < m.atom_pairs.size(); iat++) {
                    gemmi::Atom *at_1 = m.atom_pairs[iat].first;
                    gemmi::Atom *at_2 = m.atom_pairs[iat].second;
                    gemmi::Vec3 diff = at_1->pos - at_2->pos;
                    double dd = diff.length_sq();
                    sum_devi += std::sqrt(dd);
                }
                double av_devi = sum_devi / static_cast<double>(m.atom_pairs.size());
                if (av_devi < best_av_devi) {
                    best_av_devi = av_devi;
                    rm = m;
                }
            }
        }
        return rm;
    }
};

}  // namespace coot