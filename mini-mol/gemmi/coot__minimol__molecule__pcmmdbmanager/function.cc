#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>

namespace coot {
namespace minimol {

gemmi::Structure molecule::pcmmdbmanager_gemmi() const {
    gemmi::Structure structure;
    
    gemmi::Model model;
    
    for (unsigned int ifrag = 0; ifrag < fragments.size(); ifrag++) {
        const fragment& frag = fragments[ifrag];
        
        gemmi::Chain chain;
        chain.name = frag.fragment_id;
        
        for (int ires = frag.min_res_no(); ires <= frag.max_residue_number(); ires++) {
            const residue& res = frag[ires];
            if (res.atoms.size() > 0) {
                gemmi::Residue gemmi_res;
                gemmi_res.name = res.name;
                gemmi_res.seqid.num.value = res.seqnum;
                gemmi_res.seqid.icode = res.ins_code;
                
                for (unsigned int iatom = 0; iatom < res.atoms.size(); iatom++) {
                    const atom& this_atom = res.atoms[iatom];
                    
                    gemmi::Atom atom;
                    atom.pos = gemmi::Position(this_atom.pos.x, this_atom.pos.y, this_atom.pos.z);
                    atom.name = this_atom.name;
                    atom.element = gemmi::Element(this_atom.element);
                    atom.altloc = this_atom.altloc;
                    atom.occ = this_atom.occupancy;
                    atom.b_iso = this_atom.temperature_factor;
                    
                    gemmi_res.atoms.push_back(atom);
                }
                
                chain.residues.push_back(gemmi_res);
            }
        }
        
        model.chains.push_back(chain);
    }
    
    structure.models.push_back(model);
    
    if (have_cell) {
        structure.cell.a = mmdb_cell[0];
        structure.cell.b = mmdb_cell[1];
        structure.cell.c = mmdb_cell[2];
        structure.cell.alpha = mmdb_cell[3];
        structure.cell.beta = mmdb_cell[4];
        structure.cell.gamma = mmdb_cell[5];
    }
    
    if (have_spacegroup) {
        structure.spacegroup_hm = mmdb_spacegroup;
    }
    
    return structure;
}

} // namespace minimol
} // namespace coot