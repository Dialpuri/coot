#include "function.hh"

namespace coot {

std::pair<int, size_t> delete_atom_gemmi(gemmi::Structure& st, const atom_spec_t& atom_spec) {
    if (st.models.empty()) {
        return {0, 0};
    }
    
    const std::string& chain_id = atom_spec.chain_id;
    int resno = atom_spec.res_no;
    std::string ins_code = atom_spec.ins_code;
    std::string atname = atom_spec.atom_name;
    std::string altconf = atom_spec.alt_conf;
    
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_id) continue;
            
            for (auto& res : chain.residues) {
                if (res.seqid.num.value != resno) continue;
                if (res.seqid.icode != ins_code[0]) continue;
                
                for (auto it = res.atoms.begin(); it != res.atoms.end(); ++it) {
                    if (it->name != atname) continue;
                    if (it->altloc != altconf[0]) continue;
                    
                    res.atoms.erase(it);
                    
                    // Clear altloc and set occupancy to 1.0 if one atom remains
                    int n_matching = 0;
                    gemmi::Atom* remaining = nullptr;
                    for (auto& a : res.atoms) {
                        if (a.name == atname) {
                            n_matching++;
                            remaining = &a;
                        }
                    }
                    if (n_matching == 1 && remaining) {
                        remaining->altloc = '\0';
                        if (remaining->occ > 0.009) {
                            remaining->occ = 1.0;
                        }
                    }
                    
                    return {1, gemmi::count_atom_sites(st)};
                }
            }
        }
    }
    
    return {0, gemmi::count_atom_sites(st)};
}

}