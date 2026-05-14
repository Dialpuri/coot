#pragma once
#include <string>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__unused_chain_id/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__remove_TER_on_last_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__max_resno_in_chain/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__pdbcleanup_serial_residue_numbers/gemmi/function.hh"

namespace coot {
namespace molecule_t {

struct WaterAtom {
    std::string name;
    std::string element;
    gemmi::Vec3 pos;
    float temperature_factor;
    std::string altLoc;
};

struct WaterResidue {
    std::string name;
    int seqnum;
    std::vector<WaterAtom> atoms;
};

struct WaterFragment {
    std::vector<WaterResidue> residues;
};

int
insert_waters_into_molecule_gemmi(gemmi::Structure& st, const std::vector<WaterFragment>& water_mol_fragments,
                                  const std::string& res_name) {
    int water_count = 0;
    
    gemmi::Chain* solvent_chain = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "HOH" || res.name == "WAT" || res.name == " Water") {
                    solvent_chain = &chain;
                    break;
                }
            }
            if (solvent_chain) break;
        }
        if (solvent_chain) break;
    }
    
    gemmi::Chain* chain_p = solvent_chain;
    
    if (!chain_p) {
        std::pair<bool, std::string> u = coot::molecule_t::unused_chain_id_gemmi(st.models[0]);
        std::string chain_id = u.first ? u.second : "Z";
        
        gemmi::Chain new_chain;
        new_chain.name = chain_id;
        st.models[0].chains.push_back(new_chain);
        chain_p = &st.models[0].chains.back();
    } else {
        coot::molecule_t::remove_TER_on_last_residue_gemmi(*chain_p);
    }
    
    std::pair<bool, int> p = coot::util::max_resno_in_chain_gemmi(*chain_p);
    int max_resno;
    if (p.first) {
        max_resno = p.second;
    } else {
        max_resno = 0;
    }
    
    int prev_max_resno = max_resno;
    
    for (const auto& frag : water_mol_fragments) {
        for (const auto& res : frag.residues) {
            for (const auto& atom : res.atoms) {
                gemmi::Residue new_res;
                new_res.name = res_name;
                new_res.seqid.num.value = prev_max_resno + 1 + water_count;
                
                gemmi::Atom new_atom;
                new_atom.name = atom.name;
                new_atom.pos = gemmi::Position(atom.pos);
                new_atom.b_iso = atom.temperature_factor;
                new_atom.occ = 1.0;
                new_atom.element = gemmi::Element(atom.element.c_str());
                if (!atom.altLoc.empty()) {
                    new_atom.altloc = atom.altLoc[0];
                }
                
                new_res.atoms.push_back(new_atom);
                chain_p->residues.push_back(new_res);
                water_count++;
            }
        }
    }
    
    coot::util::pdbcleanup_serial_residue_numbers_gemmi(st);
    
    return 0;
}

} // namespace molecule_t
} // namespace coot