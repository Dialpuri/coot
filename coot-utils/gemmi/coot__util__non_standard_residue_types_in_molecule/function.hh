#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__residue_types_in_molecule/gemmi/function.hh"

namespace coot {
namespace util {

inline bool is_member_p(const std::vector<std::string>& v, const std::string& a) {
    bool ir = false;
    unsigned int vsize = v.size();
    for (unsigned int i = 0; i < vsize; i++) {
        if (v[i] == a) {
            ir = true;
            break;
        }
    }
    return ir;
}

inline std::vector<std::string> standard_residue_types() {
    std::vector<std::string> v;
    v.push_back("ALA"); v.push_back("ARG"); v.push_back("ASP");
    v.push_back("ASN"); v.push_back("CYS"); v.push_back("SER");
    v.push_back("PRO"); v.push_back("PHE"); v.push_back("GLY");
    v.push_back("GLU"); v.push_back("GLN"); v.push_back("ILE");
    v.push_back("LEU"); v.push_back("TYR"); v.push_back("TRP");
    v.push_back("HIS"); v.push_back("LYS"); v.push_back("MET");
    v.push_back("VAL"); v.push_back("THR"); v.push_back("MSE");
    return v;
}

inline std::vector<std::string> non_standard_residue_types_in_molecule_gemmi(const gemmi::Structure& st) {
    std::vector<std::string> r;
    std::vector<std::string> v = residue_types_in_molecule_gemmi(st);
    std::vector<std::string> standards = standard_residue_types();
    
    for (unsigned int i = 0; i < v.size(); i++) {
        if (!is_member_p(standards, v[i]))
            r.push_back(v[i]);
    }
    
    return r;
}

} // namespace util
} // namespace coot