#pragma once
#include <string>
#include <vector>
#include <optional>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// Pad a gemmi atom name to MMDB's 4-character width (PDB format)
inline std::string pad_atom_name(const std::string& name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    return name;  // already 4 chars
}

struct atom_spec_t {
    std::string chain_id;
    int res_no;
    std::string ins_code;
    std::string atom_name;
    std::string alt_conf;
    
    atom_spec_t() : chain_id(""), res_no(0), ins_code(""), atom_name(""), alt_conf("") {}
    atom_spec_t(const std::string& chain_id_, int res_no_, const std::string& ins_code_,
                const std::string& atom_name_, const std::string& alt_conf_)
        : chain_id(chain_id_), res_no(res_no_), ins_code(ins_code_),
          atom_name(atom_name_), alt_conf(alt_conf_) {}
    
    const gemmi::Atom* get_atom_gemmi(const gemmi::Structure& st) const {
        const gemmi::Atom* at = nullptr;
        if (st.models.empty()) return nullptr;
        
        for (const gemmi::Model& model : st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                if (chain.name != chain_id) continue;
                
                for (const gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value != res_no) continue;
                    
                    // Check insertion code - normalize MMDB "" to gemmi ' '
                    std::string norm_ins = ins_code.empty() ? std::string(" ") : ins_code;
                    if (std::string(1, res.seqid.icode) != norm_ins) continue;
                    
                    for (const gemmi::Atom& atom : res.atoms) {
                        // Empty alt_conf means match any alt loc
                        if (!alt_conf.empty() && atom.altloc != alt_conf[0]) continue;
                        
                        // Pad atom name to 4 chars for comparison (PDB format - right-aligned)
                        std::string padded_atom_name = pad_atom_name(atom.name);
                        if (padded_atom_name == atom_name) {
                            at = &atom;
                            return at;
                        }
                    }
                }
            }
        }
        return nullptr;
    }
};

} // namespace coot