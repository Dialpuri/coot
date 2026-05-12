#include "function.hh"
#include <stdexcept>
#include <iostream>

namespace coot {

std::string protein_geometry::get_group(const std::string& res_name_in) const {
    std::string res_name = res_name_in;
    if (res_name.empty())
        throw std::runtime_error("empty residue name");
    
    // Handle terminal caps
    if (res_name.length() > 3) {
        if (res_name.substr(0, 3) == "ACE" || res_name.substr(0, 3) == "NMA") {
            res_name = res_name.substr(3);
        }
    }
    
    // Peptide residues
    if (res_name == "ALA" || res_name == "ARG" || res_name == "ASN" ||
        res_name == "ASP" || res_name == "CYS" || res_name == "GLN" ||
        res_name == "GLU" || res_name == "GLY" || res_name == "HIS" ||
        res_name == "ILE" || res_name == "LEU" || res_name == "LYS" ||
        res_name == "MET" || res_name == "PHE" || res_name == "PRO" ||
        res_name == "SER" || res_name == "THR" || res_name == "TRP" ||
        res_name == "TYR" || res_name == "VAL") {
        return "peptide";
    }
    
    // Nucleic acids
    if (res_name == "A" || res_name == "G" || res_name == "C" ||
        res_name == "U" || res_name == "I" || res_name == "DA" ||
        res_name == "DG" || res_name == "DC" || res_name == "DT" ||
        res_name == "DI") {
        return "nucleic-acid";
    }
    
    // Water
    if (res_name == "HOH" || res_name == "WAT" || res_name == "OH2") {
        return "water";
    }
    
    // Ions
    if (res_name == "CAF" || res_name == "MG" || res_name == "CA" ||
        res_name == "ZN" || res_name == "FE" || res_name == "MN") {
        return "ion";
    }
    
    throw std::runtime_error("unknown residue type: " + res_name);
}

std::string protein_geometry::get_group_gemmi(const gemmi::Residue& residue) const {
    std::string res_name = residue.name;
    try {
        return get_group(res_name);
    }
    catch(const std::runtime_error& e) {
        std::cout << "ERROR:: pg::get_group" << res_name << " " << e.what() << '\n';
        return std::string("");
    }
}

} // namespace coot