#pragma once

#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue__832882/gemmi/function.hh"

namespace coot {
namespace util {

// Create a minimal gemmi::Structure containing just the given residue
// Takes a CRA (Chain-Residue-Atom) to preserve parent context (chain name)
// Returns a newly allocated Structure with one model, one chain containing the copied residue
inline gemmi::Structure*
create_mmdbmanager_from_residue_gemmi(const gemmi::CRA& cra) {
    if (!cra.residue || !cra.chain) {
        return nullptr;
    }
    
    gemmi::Structure* st = new gemmi::Structure;
    
    // Copy the residue with no alternate conformation filtering
    std::pair<bool, std::string> use_alt_conf(false, "");
    gemmi::Residue* rres = deep_copy_this_residue_gemmi(cra.residue, use_alt_conf);
    if (!rres) {
        delete st;
        return nullptr;
    }
    
    // Create model and chain
    gemmi::Model model;
    gemmi::Chain chain(cra.chain->name);
    chain.residues.push_back(*rres);
    model.chains.push_back(chain);
    
    st->models.push_back(model);
    
    return st;
}

} // namespace util
} // namespace coot