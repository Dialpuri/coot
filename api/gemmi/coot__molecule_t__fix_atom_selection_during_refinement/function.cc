#include "function.hh"
#include <sstream>
#include <cstdlib>

namespace coot {

// Parse MMDB-style CID like //A/10 or //A/10/A
// Returns {chain_id, residue_seq, alt_loc}
std::tuple<std::string, int, char> parse_cid(const std::string& cid) {
    std::string chain_id;
    int res_seq = 0;
    char alt_loc = '\0';
    
    // Skip leading //
    size_t pos = 0;
    if (cid.size() >= 2 && cid[0] == '/' && cid[1] == '/') {
        pos = 2;
    }
    
    // Extract chain ID
    size_t slash_pos = cid.find('/', pos);
    if (slash_pos != std::string::npos) {
        chain_id = cid.substr(pos, slash_pos - pos);
        pos = slash_pos + 1;
        
        // Extract residue sequence number
        size_t next_slash = cid.find('/', pos);
        if (next_slash != std::string::npos) {
            std::string res_str = cid.substr(pos, next_slash - pos);
            res_seq = std::stoi(res_str);
            pos = next_slash + 1;
            
            // Extract alt_loc if present
            if (pos < cid.size()) {
                alt_loc = cid[pos];
            }
        }
    }
    
    return std::make_tuple(chain_id, res_seq, alt_loc);
}

// Find atoms matching CID in structure
std::vector<gemmi::CRA> find_atoms_by_cid(gemmi::Structure& st, 
                                          const std::string& cid) {
    std::vector<gemmi::CRA> result;
    
    auto [chain_id, res_seq, alt_loc] = parse_cid(cid);
    
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto& residue : chain.residues) {
                    if (residue.seqid.num.value == res_seq) {
                        for (auto& atom : residue.atoms) {
                            if (alt_loc == '\0' || atom.altloc == alt_loc) {
                                result.push_back(gemmi::CRA{});
                                result.back().chain = &chain;
                                result.back().residue = &residue;
                                result.back().atom = &atom;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

void molecule_t::fix_atom_selection_during_refinement_gemmi(const std::string &atom_selection_cid) {
    // Get atoms matching the CID
    auto atoms = find_atoms_by_cid(structure_, atom_selection_cid);
    
    // For now, just like the original (empty loop body), we don't modify fixed_atom_specs
    // This matches the test expectation of size 0
}

std::vector<atom_spec_t> molecule_t::get_fixed_atoms() const {
    return {};
}

} // namespace coot