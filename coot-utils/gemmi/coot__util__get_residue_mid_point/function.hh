#pragma once

#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/clipper.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue__ffd79d/gemmi/function.hh"

namespace coot { namespace util {

inline std::pair<bool, clipper::Coord_orth>
get_residue_mid_point_gemmi(const std::string &chain_id, int res_no, const std::string &insertion_code, gemmi::Structure &st) {
    bool status = false;
    clipper::Coord_orth co(0,0,0);
    
    gemmi::Residue *res = coot::util::get_residue_gemmi(chain_id, res_no, insertion_code, st);
    if (!res) return std::make_pair(false, co);
    
    unsigned int n = 0;
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;
    
    for (const gemmi::Atom& atom : res->atoms) {
        sum_x += atom.pos.x;
        sum_y += atom.pos.y;
        sum_z += atom.pos.z;
        n += 1;
    }
    
    if (n > 0) {
        double nn = static_cast<double>(n);
        clipper::Coord_orth pt(sum_x/nn, sum_y/nn, sum_z/nn);
        co = pt;
        status = true;
    }
    
    return std::make_pair(status, co);
}

}} // namespace coot::util