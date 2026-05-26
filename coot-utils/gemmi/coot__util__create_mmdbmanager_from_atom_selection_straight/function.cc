#include "function.hh"
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <map>

namespace coot {
namespace util {

gemmi::Structure create_mmdbmanager_from_atom_selection_straight_gemmi(
    const gemmi::Structure& orig_st,
    const std::function<bool(const gemmi::CRA&)>& selection_callback) {
    
    gemmi::Structure result;
    result.name = orig_st.name;
    result.cell = orig_st.cell;
    result.spacegroup_hm = orig_st.spacegroup_hm;
    
    // Collect selected atoms grouped by chain, residue
    // Map: chain_name -> (seqid -> atoms)
    std::map<std::string, std::map<gemmi::SeqId, std::vector<gemmi::Atom>>> selected;
    
    for (const gemmi::Model& model : orig_st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {
                    gemmi::CRA cra{const_cast<gemmi::Chain*>(&chain), 
                                   const_cast<gemmi::Residue*>(&res), 
                                   const_cast<gemmi::Atom*>(&atom)};
                    if (selection_callback(cra)) {
                        selected[chain.name][res.seqid].push_back(atom);
                    }
                }
            }
        }
    }
    
    // Build output structure with single model
    gemmi::Model model;
    for (const auto& [chain_name, residues] : selected) {
        gemmi::Chain chain;
        chain.name = chain_name;
        for (const auto& [seqid, atoms] : residues) {
            gemmi::Residue res;
            res.seqid = seqid;
            res.atoms = atoms;
            chain.residues.push_back(res);
        }
        model.chains.push_back(chain);
    }
    result.models.push_back(model);
    
    result.connections = orig_st.connections;
    
    return result;
}

} // namespace util
} // namespace coot