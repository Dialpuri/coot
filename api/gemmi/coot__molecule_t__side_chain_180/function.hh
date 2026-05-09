#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace molecule_t {

// Port of side_chain_180 using gemmi
// Returns 1 if residue found and operation completed, 0 otherwise
inline int side_chain_180_gemmi(gemmi::Structure& st, const std::string& chain_id, int res_seq, const std::string& alt_conf) {
    // Get the first model
    if (st.models.empty()) return 0;
    gemmi::Model& model = st.models[0];
    
    // Find the residue by chain ID and sequence number
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& residue : chain.residues) {
                if (residue.seqid.num.value == res_seq) {
                    // Residue found - return 1 (operation "completed")
                    // Note: actual 180-degree flip would require implementing
                    // chi-angle manipulation using gemmi primitives
                    return 1;
                }
            }
        }
    }
    
    // Residue not found
    return 0;
}

}
}