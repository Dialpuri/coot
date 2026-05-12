#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/math.hpp>
#include <gemmi/elem.hpp>

namespace coot {
namespace ideal_rna {

gemmi::Structure ideal_rna(const std::string &RNA_or_DNA, const std::string &form,
                           short int single_stranded_flag,
                           const std::string &sequence, gemmi::Structure *standard_residues) {
    gemmi::Structure result;
    result.name = "ideal_" + RNA_or_DNA + "_" + form;
    
    gemmi::Model model;
    gemmi::Chain chain("A");
    
    std::string res_name_prefix = (RNA_or_DNA == "RNA") ? "N" : "d";
    
    for (size_t i = 0; i < sequence.length(); ++i) {
        char base = sequence[i];
        gemmi::Residue res;
        res.seqid = gemmi::SeqId(i + 1, ' ');
        
        switch (base) {
            case 'a': case 'A': res.name = res_name_prefix + "dA"; break;
            case 't': case 'T': res.name = res_name_prefix + "dT"; break;
            case 'u': case 'U': res.name = res_name_prefix + "U"; break;
            case 'g': case 'G': res.name = res_name_prefix + "dG"; break;
            case 'c': case 'C': res.name = res_name_prefix + "dC"; break;
            default: res.name = res_name_prefix + "N"; break;
        }
        
        gemmi::Atom ca; ca.name = "C1'"; ca.pos = gemmi::Position(i * 10.0, 0.0, 0.0); ca.element = gemmi::Element("C"); res.atoms.push_back(ca);
        gemmi::Atom o3; o3.name = "O3'"; o3.pos = gemmi::Position(i * 10.0 + 2.0, 1.0, 0.0); o3.element = gemmi::Element("O"); res.atoms.push_back(o3);
        gemmi::Atom p;  p.name = "P";   p.pos = gemmi::Position(i * 10.0 + 5.0, 2.0, 0.0); p.element = gemmi::Element("P"); res.atoms.push_back(p);
        
        chain.residues.push_back(res);
    }
    
    model.chains.push_back(chain);
    result.models.push_back(model);
    result.cell.set(100.0, 100.0, 100.0, 90.0, 90.0, 90.0);
    
    return result;
}

} // namespace ideal_rna
} // namespace coot