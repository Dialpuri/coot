#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>

#include "coot-utils/coot-coord-utils.hh"
#include "clipper/core/coords.h"

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <iostream>
#include <sstream>

namespace coot {
namespace util {

namespace {

static std::string trim_atom_name(const std::string& name) {
    std::string trimmed = name;
    size_t start = trimmed.find_first_not_of(' ');
    size_t end   = trimmed.find_last_not_of(' ');
    if (start == std::string::npos) return "";
    return trimmed.substr(start, end - start + 1);
}

} // anonymous namespace

inline std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth>>
get_matching_indices_gemmi(gemmi::Model* mol1,
                           gemmi::Model* mol2,
                           const coot::lsq_range_match_info_t& match_info,
                           int every_nth) {
    std::vector<clipper::Coord_orth> coords1;
    std::vector<clipper::Coord_orth> coords2;

    gemmi::Chain* chain1 = nullptr;
    for (auto& ch : mol1->chains) {
        if (ch.name == match_info.reference_chain_id) {
            chain1 = &ch;
            break;
        }
    }

    gemmi::Chain* chain2 = nullptr;
    for (auto& ch : mol2->chains) {
        if (ch.name == match_info.matcher_chain_id) {
            chain2 = &ch;
            break;
        }
    }

    if (!chain1 || !chain2) {
        return {coords1, coords2};
    }

    int n = 0;
    for (int resno = match_info.to_reference_start_resno;
         resno <= match_info.to_reference_end_resno; resno++) {

        if (n % every_nth != 0) {
            n++;
            continue;
        }
        n++;

        gemmi::Residue* res1 = nullptr;
        for (auto& r : chain1->residues) {
            if (r.seqid.num.value == resno) {
                res1 = &r;
                break;
            }
        }

        int matcher_resno = resno - match_info.to_reference_start_resno +
                            match_info.from_matcher_start_resno;
        gemmi::Residue* res2 = nullptr;
        for (auto& r : chain2->residues) {
            if (r.seqid.num.value == matcher_resno) {
                res2 = &r;
                break;
            }
        }

        if (!res1 || !res2) continue;

        if (match_info.match_type_flag == 2) {
            for (auto& a : res1->atoms) {
                std::string trimmed = trim_atom_name(a.name);
                if (trimmed == "CA") {
                    coords1.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                    break;
                }
            }
            for (auto& a : res2->atoms) {
                std::string trimmed = trim_atom_name(a.name);
                if (trimmed == "CA") {
                    coords2.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                    break;
                }
            }
        }
        else if (match_info.match_type_flag == 1) {
            for (auto& a : res1->atoms) {
                if (!a.element.is_hydrogen()) {
                    coords1.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                }
            }
            for (auto& a : res2->atoms) {
                if (!a.element.is_hydrogen()) {
                    coords2.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                }
            }
        }
        else if (match_info.match_type_flag == 3) {
            std::vector<std::string> names1;
            for (auto& a : res1->atoms) {
                if (!a.element.is_hydrogen()) {
                    coords1.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                    names1.push_back(trim_atom_name(a.name));
                }
            }
            for (auto& a : res2->atoms) {
                if (a.element.is_hydrogen()) continue;
                std::string trimmed = trim_atom_name(a.name);
                for (auto& n : names1) {
                    if (trimmed == n) {
                        coords2.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
                        break;
                    }
                }
            }
        }
    }

    return {coords1, coords2};
}

inline std::pair<short int, clipper::RTop_orth>
get_lsq_matrix_gemmi(gemmi::Model* mol1,
                     gemmi::Model* mol2,
                     const std::vector<coot::lsq_range_match_info_t>& matches,
                     int every_nth,
                     bool summary_to_screen) {

    short int istat = 0;
    clipper::RTop_orth rtop(clipper::Mat33<double>(0, 0, 0, 0, 0, 0, 0, 0, 0),
                            clipper::Coord_orth(0, 0, 0));

    std::vector<clipper::Coord_orth> co1v;
    std::vector<clipper::Coord_orth> co2v;

    for (unsigned int i = 0; i < matches.size(); i++) {
        std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth>> p =
            get_matching_indices_gemmi(mol1, mol2, matches[i], every_nth);

        if (p.first.size() > 0 && p.first.size() == p.second.size()) {
            for (unsigned int j = 0; j < p.first.size(); j++) {
                co1v.push_back(p.first[j]);
                co2v.push_back(p.second[j]);
            }
        }
    }

    if (co1v.size() > 0) {
        if (co1v.size() > 2) {
            if (co2v.size() > 2) {
                rtop = clipper::RTop_orth(co2v, co1v);
                double sum_dist  = 0.0;
                double sum_dist2 = 0.0;
                double mind  =  999999999.9;
                double maxd  = -999999999.9;
                double d;
                for (unsigned int i = 0; i < co2v.size(); i++) {
                    d = clipper::Coord_orth::length(
                        co1v[i],
                        clipper::Coord_orth(co2v[i].transform(rtop)));
                    sum_dist  += d;
                    sum_dist2 += d * d;
                    if (d > maxd) maxd = d;
                    if (d < mind) mind = d;
                }
                double mean = sum_dist / double(co2v.size());
                double v    = sum_dist2 / double(co2v.size());

                if (summary_to_screen) {
                    std::cout << "INFO:: " << co1v.size() << " matched atoms had: \n"
                              << "   mean devi: " << mean << "\n"
                              << "    rms devi: " << std::sqrt(v) << "\n"
                              << "    max devi: " << maxd << "\n"
                              << "    min devi: " << mind << std::endl;
                }
                istat = 1;
            } else {
                std::cout << "WARNING:: not enough points to do matching (matching)"
                          << std::endl;
            }
        } else {
            std::cout << "WARNING:: not enough points to do matching (reference)"
                      << std::endl;
        }
    } else {
        std::cout << "WARNING:: no points to do matching" << std::endl;
    }

    return std::pair<short int, clipper::RTop_orth>(istat, rtop);
}

} // namespace util
} // namespace coot