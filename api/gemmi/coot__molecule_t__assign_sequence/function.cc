#include "function.hh"
#include <cctype>

namespace coot { namespace molecule_t {

std::string single_letter_to_3_letter_code(char letter) {
    switch (std::toupper(static_cast<unsigned char>(letter))) {
        case 'A': return "ALA";
        case 'C': return "CYS";
        case 'D': return "ASP";
        case 'E': return "GLU";
        case 'F': return "PHE";
        case 'G': return "GLY";
        case 'H': return "HIS";
        case 'I': return "ILE";
        case 'K': return "LYS";
        case 'L': return "LEU";
        case 'M': return "MET";
        case 'N': return "ASN";
        case 'P': return "PRO";
        case 'Q': return "GLN";
        case 'R': return "ARG";
        case 'S': return "SER";
        case 'T': return "THR";
        case 'V': return "VAL";
        case 'W': return "TRP";
        case 'Y': return "TYR";
        default: return "UNK";
    }
}

bool assign_sequence_gemmi(gemmi::Model& model, const std::string& chain_id,
                          const std::string& sequence) {
    gemmi::Chain* chain = nullptr;
    for (auto& c : model.chains) {
        if (c.name == chain_id) {
            chain = &c;
            break;
        }
    }
    
    if (!chain || chain->residues.size() != sequence.size()) {
        return false;
    }
    
    for (size_t i = 0; i < chain->residues.size(); ++i) {
        chain->residues[i].name = single_letter_to_3_letter_code(sequence[i]);
    }
    
    return true;
}

}}