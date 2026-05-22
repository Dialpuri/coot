#pragma once

#include <gemmi/pdb.hpp>
#include <iostream>

namespace coot { namespace util {

// Return number of chains in the structure.
// Returns -1 if structure pointer is null or has no models.
inline int number_of_chains_gemmi(const gemmi::Structure* st) {
    if (!st)
        return -1;
    
    if (st->models.empty())
        return -1;
    
    int nchains = -1;
    for (const gemmi::Model& model : st->models) {
        nchains = model.chains.size();
        if (nchains <= 0) {
            std::cout << "bad nchains in number_of_chains_gemmi " << nchains << std::endl;
        }
    }
    return nchains;
}

}} // namespace coot::util