#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

class ShelxIns {
public:
    std::string message_for_atom_gemmi(const std::string &in_string, const gemmi::CRA &cra) const {
        std::string s = in_string;
        s += "\""; 
        s += cra.chain->name;
        s += "\""; 
        s += " "; 
        s += std::to_string(cra.residue->seqid.num.value);
        s += " ";
        s += "\""; 
        s += cra.residue->name;
        s += "\""; 
        s += " ";
        s += "\""; 
        s += cra.atom->name;
        s += "\"";
        if (cra.atom->altloc != '\0' && cra.atom->altloc != ' ') {
            s += " ,";
            s += "\""; 
            s += std::string(1, cra.atom->altloc);
            s += "\"";
        }
        return s;
    }
};

}