#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <cmath>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// Port of coot::get_environment_clash_score from MMDB to gemmi
// Calculates clash score for atoms in the environment, avoiding certain specified atoms
inline double
get_environment_clash_score_gemmi(const gemmi::Model& model,
                                  const std::vector<const gemmi::Atom*>& atoms,
                                  const std::vector<std::pair<bool, clipper::Coord_orth>>& avoid_these_atoms) {
    double cs = 0;
    double sf = 1.0;
    
    for (size_t iat = 0; iat < atoms.size(); iat++) {
        const gemmi::Atom* at = atoms[iat];
        
        // Get the residue name - need to find the parent residue
        // Since gemmi atoms don't have parent pointers, we need to search
        std::string res_name;
        bool found = false;
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& a : res.atoms) {
                    if (&a == at) {
                        res_name = res.name;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        
        // Skip ASN residues - we don't want to include clashes for ASN
        if (!found || res_name == "ASN") {
            continue;
        }
        
        clipper::Coord_orth at_pt = co_gemmi(at);
        
        for (size_t jat = 0; jat < avoid_these_atoms.size(); jat++) {
            double close_lim = 3.3;
            if (avoid_these_atoms[jat].first) close_lim = 2.5; // we can get close to waters without worry
            double close_lim_sqrd = close_lim * close_lim;
            clipper::Coord_orth avoid_pt = avoid_these_atoms[jat].second;
            double d_sqd = (at_pt - avoid_pt).lengthsq();
            
            if (d_sqd < close_lim_sqrd) {
                double diff = close_lim - std::sqrt(d_sqd);
                cs += diff * diff * sf;
            }
        }
    }
    
    return cs;
}

} // namespace coot