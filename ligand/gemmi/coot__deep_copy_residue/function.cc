#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {

gemmi::Residue deep_copy_residue_gemmi(gemmi::CRA cra) {
    gemmi::Residue res_copy;
    res_copy.name = cra.residue->name;
    res_copy.seqid.num.value = cra.residue->seqid.num.value;
    res_copy.seqid.icode = cra.residue->seqid.icode;
    
    // Copy all atoms
    for (const auto& atom : cra.residue->atoms) {
        res_copy.atoms.push_back(atom);
    }
    
    return res_copy;
}

}