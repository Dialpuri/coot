#include "function.hh"

namespace coot {
namespace molecule_t {

std::vector<double> get_residue_average_position_gemmi(
    const gemmi::Structure& st,
    const std::string& cid
) {
    std::vector<double> v;
    
    auto [chain_name, resnum] = parse_cid(cid);
    if (chain_name.empty()) {
        return v;
    }
    
    // Find the residue by chain name and sequence number
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name == chain_name) {
                for (const auto& residue : chain.residues) {
                    if (residue.seqid.num.value == resnum) {
                        // Calculate average position
                        double sum_x = 0.0;
                        double sum_y = 0.0;
                        double sum_z = 0.0;
                        int atom_count = 0;
                        
                        for (const auto& atom : residue.atoms) {
                            // Check if atom name indicates Ter (terminal)
                            // In gemmi, we can check the atom name
                            if (atom.name == "TER") {
                                continue;
                            }
                            sum_x += atom.pos.x;
                            sum_y += atom.pos.y;
                            sum_z += atom.pos.z;
                            atom_count++;
                        }
                        
                        if (atom_count > 0) {
                            double inv_count = 1.0 / static_cast<double>(atom_count);
                            v.push_back(sum_x * inv_count);
                            v.push_back(sum_y * inv_count);
                            v.push_back(sum_z * inv_count);
                        }
                        return v;
                    }
                }
            }
        }
    }
    
    // Residue not found
    return v;
}

} // namespace molecule_t
} // namespace coot